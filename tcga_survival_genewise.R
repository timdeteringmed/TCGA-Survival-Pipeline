# ============================================
# \_\_\_\_  \_  \_\_\_\_\_  \_      Rheinisch-
# \_   \_\_ \__ \_  \_  \_  \_   Westfaelische
# \_\_\_  \_\_\_\_  \_  \_\_\_      Technische
# \_  \_   \__ \__  \_  \_  \_      Hochschule
# \_   \_   \_  \_  \_  \_  \_          Aachen
# ============================================

# ==============================================================================
# Project: MD Thesis
# Focus: ITIH2 / ITIH5 Pan-Cancer Survival Analysis; REMARK-compliant Master Script
# Author: Tim Detering
# filename: TCGA_remark.R
# ==============================================================================

# ==============================================================================
# Guidelines:
#   McShane LM et al. (2005) JNCI 97:1180-1184  [REMARK]
#   Altman DG et al. (2012) PLoS Med 9:e1001216 [REMARK E&E]
#
# REMARK items addressed per section (marked inline as [REMARK #N]):
#   Item  4 – Study participants: eligibility criteria
#   Item  5 – Treatment / follow-up
#   Item  6 – Patient flow / exclusions (all reasons logged)
#   Item  7 – Assay / marker distribution reported
#   Item  8 – Cutoff pre-specified a priori, not data-driven
#   Item  9 – Number of events reported (not only N)
#   Item 10 – Multivariable model covariates pre-specified
#   Item 11 – PH assumption tested (Schoenfeld residuals)
#   Item 12 – Confidence intervals mandatory for all HRs
#   Item 13 – Missing data explicitly quantified and reported
#   Item 14 – ALL entities reported; no selective reporting
#   Item 15 – Median follow-up (reverse KM, Schemper & Smith 1996)
#   Item 16 – Univariate AND multivariate results both reported
#   Item 20 – Limitations flagged programmatically
# ==============================================================================


# --- 1. Libraries & Setup -----------------------------------------------------
suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(survival)
  library(survminer)
  library(future.apply)
  library(dplyr)
  library(readr)
  library(tibble)
})


# --- 2. Pre-specified Analysis Parameters [REMARK Item 8, 10] -----------------

TARGET_GENES  <- c("ITIH2", "ITIH5")
OUT_DIR       <- "results_tcga_REMARK"
MIN_SAMPLES   <- 30        
MIN_EVENTS    <- 10        
CUTOFF_METHOD <- "tertile"   

# [REMARK Item 10] Multivariable model covariates pre-specified:
MV_COVARIATES <- c("age_at_diagnosis", "stage_num")
MV_MIN_OBS    <- 10   

FDR_THRESHOLD <- 0.05

dir.create(OUT_DIR,                              showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "KM_plots"),       showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "exclusion_logs"), showWarnings = FALSE)


# --- 3. Helper Functions ------------------------------------------------------

# [REMARK Item 4/5] AJCC stage to ordinal integer
parse_stage <- function(x) {
  dplyr::case_when(
    grepl("Stage I(A|B|C)?$",   x, ignore.case = TRUE) ~ 1,
    grepl("Stage II(A|B|C)?$",  x, ignore.case = TRUE) ~ 2,
    grepl("Stage III(A|B|C)?$", x, ignore.case = TRUE) ~ 3,
    grepl("Stage IV(A|B|C)?$",  x, ignore.case = TRUE) ~ 4,
    TRUE ~ NA_real_
  )
}

# [REMARK Item 8] Pre-specified cutoff assignment
assign_groups <- function(expr_vec, method = CUTOFF_METHOD) {
  if (method == "median") {
    cut_val <- median(expr_vec, na.rm = TRUE)
    groups  <- ifelse(expr_vec >= cut_val, "High", "Low")
    cutoff_report <- as.character(round(cut_val, 4))
  } else if (method == "tertile") {
    q <- quantile(expr_vec, probs = c(1/3, 2/3), na.rm = TRUE)
    groups <- dplyr::case_when(
      expr_vec <= q[1] ~ "Low",
      expr_vec >= q[2] ~ "High",
      TRUE             ~ NA_character_   # middle tertile excluded [REMARK Item 6]
    )
    cutoff_report <- paste0("T1=", round(q[1], 4), "|T2=", round(q[2], 4))
  } else {
    stop("Unknown CUTOFF_METHOD: ", method)
  }
  list(groups = groups, cutoff_report = cutoff_report)
}

# [REMARK Item 15] Median follow-up via reverse Kaplan-Meier
# (Schemper & Smith 1996, Statistics in Medicine)
median_followup_reverseKM <- function(time, status) {
  tryCatch({
    fit <- survfit(Surv(time, 1 - status) ~ 1)
    stats <- surv_median(fit)
    stats$median
  }, error = function(e) NA_real_)
}

# [REMARK Item 7] Marker distribution summary
marker_distribution <- function(expr_vec, gene, project) {
  data.frame(
    Entity          = project,
    Gene            = gene,
    N_total         = length(expr_vec),
    N_zero          = sum(expr_vec == 0, na.rm = TRUE),
    Pct_zero        = round(100 * mean(expr_vec == 0, na.rm = TRUE), 1),
    Mean_log2FPKM   = round(mean(expr_vec,   na.rm = TRUE), 4),
    Median_log2FPKM = round(median(expr_vec, na.rm = TRUE), 4),
    SD_log2FPKM     = round(sd(expr_vec,     na.rm = TRUE), 4),
    Q1_log2FPKM     = round(quantile(expr_vec, 0.25, na.rm = TRUE), 4),
    Q3_log2FPKM     = round(quantile(expr_vec, 0.75, na.rm = TRUE), 4),
    stringsAsFactors = FALSE
  )
}


# --- 4. Core Analysis Function ------------------------------------------------
analyze_project <- function(project) {
  
  project_results      <- list()
  project_distributions <- list()
  
  # [REMARK Item 6] Exclusion log: track every patient removed and reason
  excl_log <- data.frame(
    Project = character(), Gene = character(),
    Stage   = character(), N_excluded = integer(),
    Reason  = character(), stringsAsFactors = FALSE
  )
  
  log_exclusion <- function(gene, stage, n_before, n_after, reason) {
    excl_log <<- rbind(excl_log, data.frame(
      Project = project, Gene = gene, Stage = stage,
      N_excluded = n_before - n_after, Reason = reason,
      stringsAsFactors = FALSE
    ))
  }
  
  message(sprintf("\n>>> [%s] Starting REMARK-compliant analysis...", project))
  
  tryCatch({
    
    # --- 4a. Data Download (project-specific dir for HPC thread-safety) -------
    gdc_dir <- file.path("GDCdata", project)
    dir.create(gdc_dir, showWarnings = FALSE, recursive = TRUE)
    
    query <- GDCquery(
      project       = project,
      data.category = "Transcriptome Profiling",
      data.type     = "Gene Expression Quantification",
      workflow.type = "STAR - Counts"
    )
    suppressMessages({
      GDCdownload(query, method = "api", directory = gdc_dir)
      se <- GDCprepare(query, directory = gdc_dir)
    })
    
    clin <- GDCquery_clinic(project = project, type = "clinical")
    rd   <- as.data.frame(rowData(se))
    n_raw_clin <- nrow(clin)
    
    # [REMARK Item 13] Quantify missing data in covariates before any merge
    missing_age   <- sum(is.na(as.numeric(clin$age_at_diagnosis)))
    missing_stage <- sum(is.na(parse_stage(clin$ajcc_pathologic_stage)))
    message(sprintf("  [%s] Missing: age=%d (%.0f%%), stage=%d (%.0f%%)",
                    project,
                    missing_age,   100 * missing_age   / nrow(clin),
                    missing_stage, 100 * missing_stage / nrow(clin)))
    
    # Clinical preparation
    clin <- clin %>%
      mutate(
        age_at_diagnosis = as.numeric(age_at_diagnosis),
        stage_num        = parse_stage(ajcc_pathologic_stage)
      )
    
    # --- 4b. Gene Loop --------------------------------------------------------
    for (gene in TARGET_GENES) {
      
      gene_idx <- which(rd$gene_name == gene)
      if (length(gene_idx) == 0) {
        message(sprintf("  [%s] SKIP %s: gene symbol not found in rowData.", project, gene))
        next
      }
      
      # Expression extraction + log2 transformation [REMARK Item 7]
      ens_id   <- rownames(rd)[gene_idx[1]]
      expr_raw <- assay(se, "fpkm_unstrand")[ens_id, ]
      
      df_gene <- data.frame(
        submitter_id = substr(names(expr_raw), 1, 12),
        expression   = log2(as.numeric(expr_raw) + 1)
      ) %>%
        group_by(submitter_id) %>%
        summarise(expression = mean(expression, na.rm = TRUE), .groups = "drop")
      
      n_expr <- nrow(df_gene)
      
      # Marker distribution report [REMARK Item 7]
      project_distributions[[gene]] <- marker_distribution(
        df_gene$expression, gene, project
      )
      
      # --- Merge with clinical ------------------------------------------------
      df_merged <- tryCatch({
        clin %>%
          inner_join(df_gene, by = "submitter_id") %>%
          mutate(
            # [REMARK Item 5] OS: primary endpoint
            OS_time   = coalesce(
              as.numeric(days_to_death),
              as.numeric(days_to_last_follow_up)
            ),
            OS_status  = as.integer(vital_status == "Dead"),
            # [REMARK Item 5] PFS: secondary endpoint
            PFS_time   = coalesce(
              as.numeric(days_to_recurrence),
              as.numeric(days_to_last_follow_up)
            ),
            PFS_status = as.integer(!is.na(days_to_recurrence))
          )
      }, error = function(e) NULL)
      
      if (is.null(df_merged)) {
        message(sprintf("  [%s] SKIP %s: merge failed.", project, gene))
        next
      }
      
      n_after_merge <- nrow(df_merged)
      log_exclusion(gene, "merge", n_raw_clin + n_expr, n_after_merge,
                    "no matching submitter_id")
      
      # [REMARK Item 6] Explicit exclusion steps with logging
      n_before <- nrow(df_merged)
      df_merged <- df_merged %>% filter(!is.na(OS_time), OS_time > 0)
      log_exclusion(gene, "OS_time_filter", n_before, nrow(df_merged),
                    "missing or zero OS_time")
      
      n_before <- nrow(df_merged)
      df_merged <- df_merged %>% filter(!is.na(expression))
      log_exclusion(gene, "expression_filter", n_before, nrow(df_merged),
                    "missing expression value")
      
      if (nrow(df_merged) < MIN_SAMPLES) {
        message(sprintf("  [%s] SKIP %s: N=%d < MIN_SAMPLES=%d after exclusions.",
                        project, gene, nrow(df_merged), MIN_SAMPLES))
        next
      }
      
      # [REMARK Item 15] Median follow-up (reverse KM)
      median_fu <- median_followup_reverseKM(df_merged$OS_time, df_merged$OS_status)
      
      # --- Group assignment [REMARK Item 8] -----------------------------------
      grp_res        <- assign_groups(df_merged$expression)
      df_merged$Group <- grp_res$groups
      
      n_before <- nrow(df_merged)
      df_merged <- df_merged %>% filter(!is.na(Group))
      log_exclusion(gene, "middle_tertile_exclusion", n_before, nrow(df_merged),
                    paste0("middle tertile excluded (", CUTOFF_METHOD, ")"))
      
      df_merged$Group <- factor(df_merged$Group, levels = c("Low", "High"))
      
      if (length(unique(df_merged$Group)) < 2) {
        message(sprintf("  [%s] SKIP %s: only one group after split.", project, gene))
        next
      }
      
      # [REMARK Item 9] Event counts
      n_events_os  <- sum(df_merged$OS_status == 1)
      n_high       <- sum(df_merged$Group == "High")
      n_low        <- sum(df_merged$Group == "Low")
      
      if (n_events_os < MIN_EVENTS) {
        message(sprintf("  [%s] SKIP %s: only %d OS events < MIN_EVENTS=%d.",
                        project, gene, n_events_os, MIN_EVENTS))
        next
      }
      
      # ---- OS Analysis -------------------------------------------------------
      fit_os       <- survfit(Surv(OS_time, OS_status) ~ Group, data = df_merged)
      p_logrank_os <- surv_pvalue(fit_os)$pval
      
      # [REMARK Item 12] Univariate Cox with CI
      cox_uv_fit  <- coxph(Surv(OS_time, OS_status) ~ Group, data = df_merged)
      cox_uv_sum  <- summary(cox_uv_fit)
      hr_uv       <- cox_uv_sum$coefficients["GroupHigh", "exp(coef)"]
      p_uv        <- cox_uv_sum$coefficients["GroupHigh", "Pr(>|z|)"]
      ci_low_uv   <- cox_uv_sum$conf.int["GroupHigh", "lower .95"]
      ci_up_uv    <- cox_uv_sum$conf.int["GroupHigh", "upper .95"]
      
      # [REMARK Item 11] Proportional hazards assumption test (Schoenfeld)
      ph_test_uv    <- tryCatch(cox.zph(cox_uv_fit),    error = function(e) NULL)
      ph_p_uv       <- if (!is.null(ph_test_uv)) ph_test_uv$table["GroupHigh", "p"] else NA_real_
      ph_violated_uv <- !is.na(ph_p_uv) && ph_p_uv < 0.05
      
      if (ph_violated_uv) {
        message(sprintf("  [REMARK Item 11 WARNING] [%s] %s OS: PH assumption violated (p=%.3f). Cox results should be interpreted with caution.",
                        project, gene, ph_p_uv))
      }
      
      # [REMARK Item 10, 16] Multivariate Cox (pre-specified covariates)
      hr_mv <- NA_real_; p_mv <- NA_real_
      ci_low_mv <- NA_real_; ci_up_mv <- NA_real_
      ph_p_mv <- NA_real_; ph_violated_mv <- NA
      
      cov_ok <- sapply(MV_COVARIATES, function(v) {
        sum(!is.na(df_merged[[v]])) >= MV_MIN_OBS
      })
      
      if (all(cov_ok)) {
        tryCatch({
          formula_mv <- as.formula(
            paste("Surv(OS_time, OS_status) ~ Group +",
                  paste(MV_COVARIATES, collapse = " + "))
          )
          cox_mv_fit  <- coxph(formula_mv, data = df_merged)
          cox_mv_sum  <- summary(cox_mv_fit)
          hr_mv       <- cox_mv_sum$coefficients["GroupHigh", "exp(coef)"]
          p_mv        <- cox_mv_sum$coefficients["GroupHigh", "Pr(>|z|)"]
          ci_low_mv   <- cox_mv_sum$conf.int["GroupHigh", "lower .95"]
          ci_up_mv    <- cox_mv_sum$conf.int["GroupHigh", "upper .95"]
          
          # [REMARK Item 11] PH test for multivariate model
          ph_test_mv    <- tryCatch(cox.zph(cox_mv_fit), error = function(e) NULL)
          ph_p_mv       <- if (!is.null(ph_test_mv)) ph_test_mv$table["GroupHigh", "p"] else NA_real_
          ph_violated_mv <- !is.na(ph_p_mv) && ph_p_mv < 0.05
          
          if (ph_violated_mv) {
            message(sprintf("  [REMARK Item 11 WARNING] [%s] %s OS MV: PH assumption violated (p=%.3f).",
                            project, gene, ph_p_mv))
          }
        }, error = function(e) {
          message(sprintf("  [%s] %s: multivariate Cox failed: %s", project, gene, e$message))
        })
      } else {
        missing_covs <- names(cov_ok)[!cov_ok]
        message(sprintf("  [REMARK Item 13] [%s] %s: MV Cox skipped – insufficient data for: %s",
                        project, gene, paste(missing_covs, collapse = ", ")))
      }
      
      # ---- PFS Analysis [REMARK Item 5] --------------------------------------
      p_logrank_pfs <- NA_real_; hr_pfs <- NA_real_; p_cox_pfs <- NA_real_
      ci_low_pfs <- NA_real_; ci_up_pfs <- NA_real_
      n_events_pfs <- NA_integer_; ph_p_pfs <- NA_real_; ph_violated_pfs <- NA
      
      df_pfs <- df_merged %>% filter(!is.na(PFS_time), PFS_time > 0)
      n_events_pfs_check <- sum(df_pfs$PFS_status == 1, na.rm = TRUE)
      
      if (nrow(df_pfs) >= MIN_SAMPLES &&
          length(unique(df_pfs$Group)) == 2 &&
          n_events_pfs_check >= MIN_EVENTS) {
        
        n_events_pfs <- n_events_pfs_check
        
        tryCatch({
          fit_pfs       <- survfit(Surv(PFS_time, PFS_status) ~ Group, data = df_pfs)
          p_logrank_pfs <- surv_pvalue(fit_pfs)$pval
          
          cox_pfs_fit  <- coxph(Surv(PFS_time, PFS_status) ~ Group, data = df_pfs)
          cox_pfs_sum  <- summary(cox_pfs_fit)
          hr_pfs       <- cox_pfs_sum$coefficients["GroupHigh", "exp(coef)"]
          p_cox_pfs    <- cox_pfs_sum$coefficients["GroupHigh", "Pr(>|z|)"]
          ci_low_pfs   <- cox_pfs_sum$conf.int["GroupHigh", "lower .95"]
          ci_up_pfs    <- cox_pfs_sum$conf.int["GroupHigh", "upper .95"]
          
          # [REMARK Item 11] PH test for PFS
          ph_test_pfs    <- tryCatch(cox.zph(cox_pfs_fit), error = function(e) NULL)
          ph_p_pfs       <- if (!is.null(ph_test_pfs)) ph_test_pfs$table["GroupHigh", "p"] else NA_real_
          ph_violated_pfs <- !is.na(ph_p_pfs) && ph_p_pfs < 0.05
          
          # PFS KM plot
          pdf(file.path(OUT_DIR, "KM_plots",
                        paste0(project, "_", gene, "_PFS_KM.pdf")),
              width = 8, height = 7, onefile = FALSE)
          print(ggsurvplot(fit_pfs, data = df_pfs,
                           pval = TRUE, risk.table = TRUE,
                           conf.int = TRUE,
                           title       = paste(project, "-", gene, "| PFS"),
                           palette     = c("#2E9FDF", "#E7B800"),
                           legend.labs = c("Low", "High")))
          dev.off()
        }, error = function(e) {
          message(sprintf("  [%s] %s: PFS analysis failed: %s", project, gene, e$message))
        })
      }
      
      # ---- OS KM plot [REMARK Items 9, 12] -----------------------------------
      pdf(file.path(OUT_DIR, "KM_plots",
                    paste0(project, "_", gene, "_OS_KM.pdf")),
          width = 8, height = 7, onefile = FALSE)
      print(ggsurvplot(fit_os, data = df_merged,
                       pval = TRUE, risk.table = TRUE,
                       conf.int = TRUE,       # [REMARK Item 12]
                       title       = paste(project, "-", gene, "| OS"),
                       palette     = c("#2E9FDF", "#E7B800"),
                       legend.labs = c("Low", "High")))
      dev.off()
      
      # ---- [REMARK Item 20] Limitations flag ---------------------------------
      limitations <- c()
      if (ph_violated_uv)   limitations <- c(limitations, "PH_violated_univariate")
      if (isTRUE(ph_violated_mv))  limitations <- c(limitations, "PH_violated_multivariate")
      if (isTRUE(ph_violated_pfs)) limitations <- c(limitations, "PH_violated_PFS")
      if (is.na(hr_mv))     limitations <- c(limitations, "multivariate_not_estimable")
      if (missing_age / n_raw_clin > 0.20)
        limitations <- c(limitations, "age_missing_>20pct")
      if (missing_stage / n_raw_clin > 0.20)
        limitations <- c(limitations, "stage_missing_>20pct")
      limitations_str <- if (length(limitations) > 0) paste(limitations, collapse = "; ") else "none"
      
      # ---- Collecting results [REMARK Items 9, 12, 14, 15, 16] -----------------
      project_results[[gene]] <- data.frame(
        # Identity
        Entity                 = project,
        Gene                   = gene,
        Cutoff_method          = CUTOFF_METHOD,
        Cutoff_value           = grp_res$cutoff_report,
        # [REMARK Item 6] Evaluable patients
        N_evaluable            = nrow(df_merged),
        N_High                 = n_high,
        N_Low                  = n_low,
        # [REMARK Item 9] Event counts
        N_events_OS            = n_events_os,
        N_events_PFS           = ifelse(is.na(n_events_pfs), NA_integer_, n_events_pfs),
        # [REMARK Item 15] Follow-up
        Median_followup_days   = round(median_fu, 1),
        # [REMARK Items 12, 16] OS univariate
        HR_OS_uv               = round(hr_uv,     4),
        CI95_low_OS_uv         = round(ci_low_uv, 4),
        CI95_up_OS_uv          = round(ci_up_uv,  4),
        P_OS_univariate        = p_uv,
        P_OS_LogRank           = p_logrank_os,
        # [REMARK Item 11] PH test OS univariate
        PH_p_OS_uv             = round(ph_p_uv, 4),
        PH_violated_OS_uv      = ph_violated_uv,
        # [REMARK Items 10, 12, 16] OS multivariate
        HR_OS_mv               = round(hr_mv,     4),
        CI95_low_OS_mv         = round(ci_low_mv, 4),
        CI95_up_OS_mv          = round(ci_up_mv,  4),
        P_OS_multivariate      = p_mv,
        PH_p_OS_mv             = round(ph_p_mv, 4),
        PH_violated_OS_mv      = ph_violated_mv,
        # [REMARK Item 5] PFS
        HR_PFS                 = round(hr_pfs,     4),
        CI95_low_PFS           = round(ci_low_pfs, 4),
        CI95_up_PFS            = round(ci_up_pfs,  4),
        P_PFS_LogRank          = p_logrank_pfs,
        P_PFS_Cox              = p_cox_pfs,
        PH_p_PFS               = round(ph_p_pfs, 4),
        PH_violated_PFS        = ph_violated_pfs,
        # [REMARK Item 20] Limitations
        Limitations            = limitations_str,
        stringsAsFactors       = FALSE
      )
      
    } # end gene loop
    
    # Write exclusion log [REMARK Item 6]
    if (nrow(excl_log) > 0) {
      write_csv(excl_log,
                file.path(OUT_DIR, "exclusion_logs",
                          paste0(project, "_exclusion_log.csv")))
    }
    
    rm(se, clin, rd)
    gc()
    
    return(list(
      results       = bind_rows(project_results),
      distributions = bind_rows(project_distributions)
    ))
    
  }, error = function(e) {
    message(sprintf("[%s] FATAL ERROR: %s", project, e$message))
    return(NULL)
  })
}


# --- 5. Execute ---------------------------------------------------------------
n_cores <- suppressWarnings(as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK")))
if (is.na(n_cores) || n_cores < 1) n_cores <- 4
plan(multicore, workers = n_cores)
message(sprintf("Using %d cores.", n_cores))

all_projects   <- grep("^TCGA-", getGDCprojects()$project_id, value = TRUE)
message(sprintf("Submitting %d TCGA projects...", length(all_projects)))

raw_results <- future_lapply(all_projects, analyze_project, future.seed = TRUE)


# --- 6. Post-processing -------------------------------------------------------

# [REMARK Item 14] Collect ALL entities (including non-significant)
final_df <- bind_rows(lapply(raw_results, function(x) if (!is.null(x)) x$results))
dist_df  <- bind_rows(lapply(raw_results, function(x) if (!is.null(x)) x$distributions))

# [REMARK Item 14] FDR correction per gene across all tested entities
final_df <- final_df %>%
  group_by(Gene) %>%
  mutate(
    FDR_OS_LogRank      = p.adjust(P_OS_LogRank,      method = "BH"),
    FDR_OS_univariate   = p.adjust(P_OS_univariate,   method = "BH"),
    FDR_OS_multivariate = p.adjust(P_OS_multivariate, method = "BH"),
    FDR_PFS_LogRank     = p.adjust(P_PFS_LogRank,     method = "BH"),
    Sig_OS_uv           = !is.na(FDR_OS_LogRank)    & FDR_OS_LogRank    < FDR_THRESHOLD,
    Sig_OS_mv           = !is.na(FDR_OS_multivariate) & FDR_OS_multivariate < FDR_THRESHOLD,
    Sig_PFS             = !is.na(FDR_PFS_LogRank)   & FDR_PFS_LogRank   < FDR_THRESHOLD,
    Direction_OS = dplyr::case_when(
      Sig_OS_uv & HR_OS_uv > 1 ~ "poor_prognosis",
      Sig_OS_uv & HR_OS_uv < 1 ~ "good_prognosis",
      TRUE                     ~ "ns"
    )
  ) %>%
  ungroup() %>%
  arrange(Gene, FDR_OS_LogRank)


# --- 7. Export ----------------------------------------------------------------

# [REMARK Item 14] Primary output: ALL entities, no filtering
write_csv(final_df,
          file.path(OUT_DIR, "REMARK_PanCancer_Summary_ALL.csv"))

# Marker distributions [REMARK Item 7]
write_csv(dist_df,
          file.path(OUT_DIR, "REMARK_Marker_Distributions.csv"))

# Convenience subset: significant hits (for downstream analyses only)
# NOTE [REMARK Item 14]: This file must NOT be the sole reported output.
sig_df <- final_df %>%
  filter(Sig_OS_uv | Sig_OS_mv | Sig_PFS)
write_csv(sig_df,
          file.path(OUT_DIR, "REMARK_Significant_Hits_FDR05.csv"))

# Entities with PH violations [REMARK Item 11]
ph_flags <- final_df %>%
  filter(isTRUE(PH_violated_OS_uv) | isTRUE(PH_violated_OS_mv) | isTRUE(PH_violated_PFS))
if (nrow(ph_flags) > 0) {
  write_csv(ph_flags, file.path(OUT_DIR, "REMARK_PH_Violations.csv"))
  message(sprintf("\n[REMARK Item 11] %d entity-gene pairs with PH violations written to REMARK_PH_Violations.csv",
                  nrow(ph_flags)))
}

# Consolidated exclusion log [REMARK Item 6]
excl_files <- list.files(file.path(OUT_DIR, "exclusion_logs"), full.names = TRUE)
if (length(excl_files) > 0) {
  all_excl <- bind_rows(lapply(excl_files, read_csv, show_col_types = FALSE))
  write_csv(all_excl, file.path(OUT_DIR, "REMARK_Exclusion_Log_ALL.csv"))
}

message(sprintf(
  "\n=== REMARK Analysis complete ===
  Total entity-gene pairs evaluated : %d
  Significant OS (FDR<%.2f, log-rank): %d
  Significant OS (FDR<%.2f, MV Cox)  : %d
  Significant PFS (FDR<%.2f)         : %d
  PH assumption violations           : %d
  Output directory                   : %s",
  nrow(final_df),
  FDR_THRESHOLD, sum(final_df$Sig_OS_uv,  na.rm = TRUE),
  FDR_THRESHOLD, sum(final_df$Sig_OS_mv,  na.rm = TRUE),
  FDR_THRESHOLD, sum(final_df$Sig_PFS,    na.rm = TRUE),
  nrow(ph_flags),
  OUT_DIR
))
