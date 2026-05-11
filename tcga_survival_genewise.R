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
# Outputs:
#   - tables/MASTER_PanCancer_Survival_REMARK.csv
#   - tables/Marker_Distributions.csv
#   - tables/Exclusion_Log_consolidated.csv
#   - KM_plots/<entity>_<gene>_OS_KM.pdf, _PFS_KM.pdf
#   - raw_plot_data/<entity>_<gene>_raw.csv
#
# Guidelines:
#   McShane LM et al. (2005) JNCI 97:1180-1184  [REMARK]
#   Altman DG et al. (2012) PLoS Med 9:e1001216 [REMARK E&E]
# ==============================================================================


# --- 1. Libraries -------------------------------------------------------------
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


# --- 2. HPC-Paths -------------------------------------------------------------

HPC_USER       <- Sys.getenv("USER", unset = "userid")
HPCWORK_BASE   <- file.path("/hpcwork", HPC_USER)
WORK_BASE      <- file.path("/work",    HPC_USER)

DOWNLOAD_BASE <- if (dir.exists(HPCWORK_BASE) && file.access(HPCWORK_BASE, mode = 2) == 0) {
  file.path(HPCWORK_BASE, "tcga_remark", "GDCdata")
} else {
  message("WARNUNG: /hpcwork nicht verfuegbar - falle zurueck auf ./GDCdata")
  "GDCdata"
}

OUT_DIR <- if (dir.exists(WORK_BASE) && file.access(WORK_BASE, mode = 2) == 0) {
  file.path(WORK_BASE, "results_tcga_REMARK")
} else {
  message("WARNUNG: /work nicht verfuegbar - falle zurueck auf ./results_tcga_REMARK")
  "results_tcga_REMARK"
}

dir.create(DOWNLOAD_BASE,                         showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_DIR,                               showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "KM_plots"),        showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "exclusion_logs"),  showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "raw_plot_data"),   showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "tables"),          showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "diagnostics"),     showWarnings = FALSE)

message(sprintf("DOWNLOAD_BASE: %s", DOWNLOAD_BASE))
message(sprintf("OUT_DIR:       %s", OUT_DIR))


# --- 3. Pre-specified Analysis Parameters [REMARK Item 8, 10] -----------------
TARGET_GENES      <- c("ITIH2", "ITIH5")
MIN_SAMPLES       <- 30
MIN_EVENTS        <- 10
CUTOFF_METHOD     <- "tertile"
MV_COVARIATES     <- c("age_at_diagnosis", "stage_num")
MV_MIN_OBS        <- 10
FDR_THRESHOLD     <- 0.05
CLEANUP_AFTER_RUN <- FALSE


# --- 4. Helper Functions ------------------------------------------------------

parse_stage <- function(x) {
  if (is.null(x) || length(x) == 0) return(numeric(0))
  x_chr <- as.character(x)
  out <- dplyr::case_when(
    grepl("Stage I(A|B|C)?$",   x_chr, ignore.case = TRUE) ~ 1.0,
    grepl("Stage II(A|B|C)?$",  x_chr, ignore.case = TRUE) ~ 2.0,
    grepl("Stage III(A|B|C)?$", x_chr, ignore.case = TRUE) ~ 3.0,
    grepl("Stage IV(A|B|C)?$",  x_chr, ignore.case = TRUE) ~ 4.0,
    TRUE ~ NA_real_
  )
  as.numeric(out)
}

assign_groups <- function(expr_vec, method = CUTOFF_METHOD) {
  if (method == "median") {
    cut_val <- median(expr_vec, na.rm = TRUE)
    return(list(groups = ifelse(expr_vec >= cut_val, "High", "Low"),
                cutoff_report = as.character(round(cut_val, 4))))
  }
  if (method == "tertile") {
    q <- quantile(expr_vec, probs = c(1/3, 2/3), na.rm = TRUE)
    grps <- dplyr::case_when(
      expr_vec <= q[1] ~ "Low",
      expr_vec >= q[2] ~ "High",
      TRUE             ~ NA_character_
    )
    return(list(groups = grps,
                cutoff_report = paste0("T1=", round(q[1], 4), "|T2=", round(q[2], 4))))
  }
  stop("Unknown CUTOFF_METHOD: ", method)
}

median_followup_reverseKM <- function(time, status) {
  tryCatch({
    fit <- survfit(Surv(time, 1 - status) ~ 1)
    surv_median(fit)$median
  }, error = function(e) NA_real_)
}

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

safe_row_value <- function(mat, row_pattern, col_name, default = NA_real_) {
  if (is.null(mat)) return(default)
  rn <- rownames(mat)
  cn <- colnames(mat)
  if (is.null(rn) || is.null(cn)) return(default)
  if (!col_name %in% cn) return(default)

  row_idx <- which(rn == row_pattern)
  if (length(row_idx) == 0) {
    # Fallback: per grep
    row_idx <- grep(row_pattern, rn, fixed = TRUE)
  }
  if (length(row_idx) == 0) return(default)
  as.numeric(mat[row_idx[1], col_name])
}


# --- 5. Phase 1 ---------------------------------------

download_project <- function(project) {
  gdc_dir <- file.path(DOWNLOAD_BASE, project)
  dir.create(gdc_dir, showWarnings = FALSE, recursive = TRUE)
  message(sprintf("[%s] Download starting...", project))

  query <- tryCatch({
    GDCquery(project = project,
             data.category = "Transcriptome Profiling",
             data.type     = "Gene Expression Quantification",
             workflow.type = "STAR - Counts")
  }, error = function(e) {
    message(sprintf("[%s] GDCquery failed: %s", project, e$message))
    NULL
  })
  if (is.null(query)) return(FALSE)

  for (attempt in 1:3) {
    ok <- tryCatch({
      suppressMessages(GDCdownload(query, method = "api", directory = gdc_dir))
      TRUE
    }, error = function(e) {
      message(sprintf("[%s] Download attempt %d/3 failed: %s",
                      project, attempt, e$message))
      if (attempt < 3) {
        unlink(gdc_dir, recursive = TRUE)
        dir.create(gdc_dir, showWarnings = FALSE, recursive = TRUE)
        Sys.sleep(15)
      }
      FALSE
    })
    if (ok) {
      message(sprintf("[%s] Download OK.", project))
      return(TRUE)
    }
  }
  return(FALSE)
}


# --- 6. Phase 2 ----------------------------------------------

analyze_gene <- function(project, gene, se, clin, rd, n_raw_clin) {

  excl_log_local <- data.frame()
  add_excl <- function(stage, n_before, n_after, reason) {
    excl_log_local <<- rbind(excl_log_local, data.frame(
      Project = project, Gene = gene, Stage = stage,
      N_excluded = n_before - n_after, Reason = reason,
      stringsAsFactors = FALSE
    ))
  }

  step_name    <- "init"
  step_context <- list()   

  set_step <- function(name, ctx = list()) {
    step_name    <<- name
    step_context <<- ctx
  }

  out <- tryCatch({

    # --- Gene Lookup --------------------------------------------------------
    set_step("gene_lookup", list(target_gene = gene, n_genes_in_rd = nrow(rd)))
    gene_idx <- which(rd$gene_name == gene)
    if (length(gene_idx) == 0) {
      message(sprintf("  [%s] SKIP %s: gene not found in rowData.", project, gene))
      return(list(result = NULL, dist = NULL, excl = excl_log_local))
    }

    ens_id <- rownames(rd)[gene_idx[1]]

    # --- Assay-Testing -----------------------------------------------------
    set_step("assay_check", list(assay_names = paste(assayNames(se), collapse = ",")))
    if (!"fpkm_unstrand" %in% assayNames(se)) {
      available <- paste(assayNames(se), collapse = ", ")
      message(sprintf("  [%s] SKIP %s: assay 'fpkm_unstrand' missing. Available: %s",
                      project, gene, available))
      return(list(result = NULL, dist = NULL, excl = excl_log_local))
    }

    # --- Expression-Extraction ----------------------------------------------
    set_step("expression_extraction", list(ens_id = ens_id))
    expr_mat <- assay(se, "fpkm_unstrand")
    if (!ens_id %in% rownames(expr_mat)) {
      message(sprintf("  [%s] SKIP %s: ens_id %s not in assay rownames.",
                      project, gene, ens_id))
      return(list(result = NULL, dist = NULL, excl = excl_log_local))
    }
    expr_raw <- expr_mat[ens_id, ]

    set_step("df_gene_creation", list(n_samples = length(expr_raw)))
    df_gene <- data.frame(
      submitter_id = substr(names(expr_raw), 1, 12),
      expression   = log2(as.numeric(expr_raw) + 1),
      stringsAsFactors = FALSE
    ) %>%
      group_by(submitter_id) %>%
      summarise(expression = mean(expression, na.rm = TRUE), .groups = "drop")

    set_step("marker_distribution", list(n_unique_patients = nrow(df_gene)))
    dist_row <- marker_distribution(df_gene$expression, gene, project)

    # --- Merge --------------------------------------------------------------
    set_step("clinical_merge",
             list(n_clin = nrow(clin), n_gene = nrow(df_gene),
                  clin_cols = paste(head(colnames(clin), 5), collapse = ",")))
    merged <- clin %>% inner_join(df_gene, by = "submitter_id")

    set_step("required_columns_check",
             list(merged_cols = paste(head(colnames(merged), 10), collapse = ",")))
    required_cols <- c("days_to_last_follow_up", "vital_status")
    missing_req   <- setdiff(required_cols, colnames(merged))
    if (length(missing_req) > 0) {
      message(sprintf("  [%s] SKIP %s: missing required cols: %s",
                      project, gene, paste(missing_req, collapse = ", ")))
      return(list(result = NULL, dist = dist_row, excl = excl_log_local))
    }

    set_step("endpoint_construction",
             list(has_recurrence = "days_to_recurrence" %in% colnames(merged),
                  has_death      = "days_to_death"      %in% colnames(merged),
                  n_merged       = nrow(merged)))
    has_recurrence <- "days_to_recurrence" %in% colnames(merged)
    has_death      <- "days_to_death"      %in% colnames(merged)

    days_death <- if (has_death) suppressWarnings(as.numeric(merged$days_to_death)) else rep(NA_real_, nrow(merged))
    days_fu    <- suppressWarnings(as.numeric(merged$days_to_last_follow_up))
    days_rec   <- if (has_recurrence) suppressWarnings(as.numeric(merged$days_to_recurrence)) else rep(NA_real_, nrow(merged))

    merged$OS_time    <- ifelse(!is.na(days_death), days_death, days_fu)
    merged$OS_status  <- as.integer(merged$vital_status == "Dead")
    merged$PFS_time   <- ifelse(!is.na(days_rec),   days_rec,   days_fu)
    merged$PFS_status <- if (has_recurrence) as.integer(!is.na(days_rec)) else 0L

    df_merged <- merged

    # --- Filter w/ Logging -------------------------------------------------
    set_step("filter_OS_time", list(n_before = nrow(df_merged)))
    n0 <- nrow(df_merged)
    df_merged <- df_merged[!is.na(df_merged$OS_time) & df_merged$OS_time > 0, ]
    add_excl("OS_time_filter", n0, nrow(df_merged), "OS_time NA or <= 0")

    set_step("filter_expression", list(n_before = nrow(df_merged)))
    n0 <- nrow(df_merged)
    df_merged <- df_merged[!is.na(df_merged$expression), ]
    add_excl("expression_filter", n0, nrow(df_merged), "expression NA")

    set_step("min_samples_check", list(n_after_filters = nrow(df_merged), min_required = MIN_SAMPLES))
    if (nrow(df_merged) < MIN_SAMPLES) {
      message(sprintf("  [%s] SKIP %s: N=%d < %d.",
                      project, gene, nrow(df_merged), MIN_SAMPLES))
      return(list(result = NULL, dist = dist_row, excl = excl_log_local))
    }

    set_step("median_followup", list(n = nrow(df_merged),
                                      n_events = sum(df_merged$OS_status == 1, na.rm = TRUE)))
    median_fu <- median_followup_reverseKM(df_merged$OS_time, df_merged$OS_status)

    # --- Cutoff ---------------------------------------------------
    set_step("group_assignment", list(method = CUTOFF_METHOD,
                                       expr_min = min(df_merged$expression, na.rm = TRUE),
                                       expr_max = max(df_merged$expression, na.rm = TRUE)))
    grp_res         <- assign_groups(df_merged$expression)
    df_merged$Group <- grp_res$groups

    set_step("filter_middle_tertile", list(n_before = nrow(df_merged)))
    n0 <- nrow(df_merged)
    df_merged <- df_merged[!is.na(df_merged$Group), ]
    add_excl("tertile_middle", n0, nrow(df_merged),
             sprintf("middle tertile excluded (%s)", CUTOFF_METHOD))

    set_step("factor_creation",
             list(unique_groups = paste(unique(df_merged$Group), collapse = ",")))
    df_merged$Group <- factor(df_merged$Group, levels = c("Low", "High"))

    if (length(unique(na.omit(df_merged$Group))) < 2) {
      message(sprintf("  [%s] SKIP %s: single group only.", project, gene))
      return(list(result = NULL, dist = dist_row, excl = excl_log_local))
    }

    set_step("events_check",
             list(n_events_os = sum(df_merged$OS_status == 1, na.rm = TRUE),
                  n_high      = sum(df_merged$Group == "High", na.rm = TRUE),
                  n_low       = sum(df_merged$Group == "Low",  na.rm = TRUE)))
    n_events_os <- sum(df_merged$OS_status == 1, na.rm = TRUE)
    if (n_events_os < MIN_EVENTS) {
      message(sprintf("  [%s] SKIP %s: only %d OS events.", project, gene, n_events_os))
      return(list(result = NULL, dist = dist_row, excl = excl_log_local))
    }

    # --- OS Univariate ------------------------------------------------------
    set_step("OS_survfit", list(n = nrow(df_merged), n_events = n_events_os))
    fit_os <- survfit(Surv(OS_time, OS_status) ~ Group, data = df_merged)

    set_step("OS_logrank")
    p_logrank_os <- surv_pvalue(fit_os, data = df_merged)$pval

    set_step("OS_cox_univariate")
    cox_uv_fit <- coxph(Surv(OS_time, OS_status) ~ Group, data = df_merged)
    cox_uv     <- summary(cox_uv_fit)

    set_step("OS_cox_uv_extract",
             list(coef_rownames = paste(rownames(cox_uv$coefficients), collapse = ","),
                  conf_rownames = paste(rownames(cox_uv$conf.int), collapse = ",")))
    hr_uv      <- safe_row_value(cox_uv$coefficients, "GroupHigh", "exp(coef)")
    p_uv       <- safe_row_value(cox_uv$coefficients, "GroupHigh", "Pr(>|z|)")
    ci_low_uv  <- safe_row_value(cox_uv$conf.int,     "GroupHigh", "lower .95")
    ci_up_uv   <- safe_row_value(cox_uv$conf.int,     "GroupHigh", "upper .95")

    set_step("OS_cox_zph_uv")
    ph_test_uv <- tryCatch(cox.zph(cox_uv_fit), error = function(e) NULL)
    ph_p_uv    <- if (!is.null(ph_test_uv))
                     safe_row_value(ph_test_uv$table, "GroupHigh", "p")
                  else NA_real_
    ph_violated_uv <- !is.na(ph_p_uv) && ph_p_uv < 0.05

    # --- OS Multivariate ----------------------------------------------------
    set_step("OS_cox_multivariate_check",
             list(covariates = paste(MV_COVARIATES, collapse = ",")))
    hr_mv <- ci_low_mv <- ci_up_mv <- p_mv <- ph_p_mv <- NA_real_
    ph_violated_mv <- NA

    cov_ok <- vapply(MV_COVARIATES,
                     function(v) v %in% colnames(df_merged) &&
                                 sum(!is.na(df_merged[[v]])) >= MV_MIN_OBS,
                     logical(1))

    if (all(cov_ok)) {
      set_step("OS_cox_multivariate_fit")
      tryCatch({
        fmla <- as.formula(paste("Surv(OS_time, OS_status) ~ Group +",
                                  paste(MV_COVARIATES, collapse = " + ")))
        cox_mv <- coxph(fmla, data = df_merged)
        smv <- summary(cox_mv)
        hr_mv     <- safe_row_value(smv$coefficients, "GroupHigh", "exp(coef)")
        p_mv      <- safe_row_value(smv$coefficients, "GroupHigh", "Pr(>|z|)")
        ci_low_mv <- safe_row_value(smv$conf.int,     "GroupHigh", "lower .95")
        ci_up_mv  <- safe_row_value(smv$conf.int,     "GroupHigh", "upper .95")

        ph_mv <- tryCatch(cox.zph(cox_mv), error = function(e) NULL)
        if (!is.null(ph_mv)) {
          ph_p_mv        <- safe_row_value(ph_mv$table, "GroupHigh", "p")
          ph_violated_mv <- !is.na(ph_p_mv) && ph_p_mv < 0.05
        }
      }, error = function(e) {
        message(sprintf("  [%s] %s: MV Cox failed: %s", project, gene, e$message))
      })
    }

    # --- PFS ----------------------------------------------------------------
    set_step("PFS_setup")
    hr_pfs <- p_logrank_pfs <- p_cox_pfs <- ci_low_pfs <- ci_up_pfs <- NA_real_
    ph_p_pfs        <- NA_real_
    ph_violated_pfs <- NA
    n_events_pfs    <- NA_integer_

    df_pfs <- df_merged[!is.na(df_merged$PFS_time) & df_merged$PFS_time > 0, ]
    n_ev_pfs_check <- sum(df_pfs$PFS_status == 1, na.rm = TRUE)

    if (nrow(df_pfs) >= MIN_SAMPLES &&
        length(unique(na.omit(df_pfs$Group))) == 2 &&
        n_ev_pfs_check >= MIN_EVENTS) {

      n_events_pfs <- n_ev_pfs_check
      set_step("PFS_analysis", list(n_pfs = nrow(df_pfs), n_events_pfs = n_ev_pfs_check))
      tryCatch({
        fit_pfs       <- survfit(Surv(PFS_time, PFS_status) ~ Group, data = df_pfs)
        p_logrank_pfs <- surv_pvalue(fit_pfs, data = df_pfs)$pval

        cox_pfs   <- coxph(Surv(PFS_time, PFS_status) ~ Group, data = df_pfs)
        spfs      <- summary(cox_pfs)
        hr_pfs    <- safe_row_value(spfs$coefficients, "GroupHigh", "exp(coef)")
        p_cox_pfs <- safe_row_value(spfs$coefficients, "GroupHigh", "Pr(>|z|)")
        ci_low_pfs<- safe_row_value(spfs$conf.int,     "GroupHigh", "lower .95")
        ci_up_pfs <- safe_row_value(spfs$conf.int,     "GroupHigh", "upper .95")

        ph_pfs <- tryCatch(cox.zph(cox_pfs), error = function(e) NULL)
        if (!is.null(ph_pfs)) {
          ph_p_pfs        <- safe_row_value(ph_pfs$table, "GroupHigh", "p")
          ph_violated_pfs <- !is.na(ph_p_pfs) && ph_p_pfs < 0.05
        }

        set_step("PFS_KM_plot")
        tryCatch({
          pdf(file.path(OUT_DIR, "KM_plots",
                        sprintf("%s_%s_PFS_KM.pdf", project, gene)),
              width = 8, height = 7, onefile = FALSE)
          print(ggsurvplot(fit_pfs, data = df_pfs,
                           pval = TRUE, risk.table = TRUE, conf.int = TRUE,
                           title       = sprintf("%s - %s | PFS", project, gene),
                           palette     = c("#2E9FDF", "#E7B800"),
                           legend.labs = c("Low", "High")))
          dev.off()
        }, error = function(e) {
          tryCatch(dev.off(), error = function(e2) NULL)
          message(sprintf("  [%s] %s: PFS plot failed: %s", project, gene, e$message))
        })
      }, error = function(e) {
        message(sprintf("  [%s] %s: PFS analysis failed: %s",
                        project, gene, e$message))
      })
    }

    # --- OS KM-Plot ---------------------------------------------------------
    set_step("OS_KM_plot")
    tryCatch({
      pdf(file.path(OUT_DIR, "KM_plots",
                    sprintf("%s_%s_OS_KM.pdf", project, gene)),
          width = 8, height = 7, onefile = FALSE)
      print(ggsurvplot(fit_os, data = df_merged,
                       pval = TRUE, risk.table = TRUE, conf.int = TRUE,
                       title       = sprintf("%s - %s | OS", project, gene),
                       palette     = c("#2E9FDF", "#E7B800"),
                       legend.labs = c("Low", "High")))
      dev.off()
    }, error = function(e) {
      tryCatch(dev.off(), error = function(e2) NULL)
      message(sprintf("  [%s] %s: OS plot failed: %s", project, gene, e$message))
    })

    # --- Raw-Data CSV -------------------------------------------------------
    set_step("raw_data_csv_write")
    tryCatch({
      keep_cols <- intersect(c("submitter_id", "expression", "Group",
                                "OS_time", "OS_status", "PFS_time", "PFS_status",
                                "age_at_diagnosis", "stage_num"),
                              colnames(df_merged))
      raw_data <- df_merged[, keep_cols, drop = FALSE]
      raw_data$Entity <- project
      raw_data$Gene   <- gene
      write_csv(raw_data,
                file.path(OUT_DIR, "raw_plot_data",
                          sprintf("%s_%s_raw.csv", project, gene)))
    }, error = function(e) {
      message(sprintf("  [%s] %s: raw CSV write failed: %s", project, gene, e$message))
    })

    set_step("results_collection")

    # --- Limitations [REMARK Item 20] --------------------------------------
    limitations <- c()
    if (ph_violated_uv)            limitations <- c(limitations, "PH_violated_OS_uv")
    if (isTRUE(ph_violated_mv))    limitations <- c(limitations, "PH_violated_OS_mv")
    if (isTRUE(ph_violated_pfs))   limitations <- c(limitations, "PH_violated_PFS")
    if (is.na(hr_mv))              limitations <- c(limitations, "MV_not_estimable")
    lim_str <- if (length(limitations) > 0) paste(limitations, collapse = ";") else "none"

    # --- Result-Rows ------------------------------------------------
    result_row <- data.frame(
      Entity                = project,
      Gene                  = gene,
      Cutoff_method         = CUTOFF_METHOD,
      Cutoff_value          = grp_res$cutoff_report,
      N_evaluable           = nrow(df_merged),
      N_High                = sum(df_merged$Group == "High", na.rm = TRUE),
      N_Low                 = sum(df_merged$Group == "Low",  na.rm = TRUE),
      N_events_OS           = n_events_os,
      N_events_PFS          = ifelse(is.na(n_events_pfs), NA_integer_, n_events_pfs),
      Median_followup_days  = round(median_fu, 1),
      HR_OS_uv              = round(hr_uv,    4),
      CI95_low_OS_uv        = round(ci_low_uv, 4),
      CI95_up_OS_uv         = round(ci_up_uv,  4),
      P_OS_univariate       = p_uv,
      P_OS_LogRank          = p_logrank_os,
      PH_p_OS_uv            = round(ph_p_uv, 4),
      PH_violated_OS_uv     = ph_violated_uv,
      HR_OS_mv              = round(hr_mv,    4),
      CI95_low_OS_mv        = round(ci_low_mv, 4),
      CI95_up_OS_mv         = round(ci_up_mv,  4),
      P_OS_multivariate     = p_mv,
      PH_p_OS_mv            = round(ph_p_mv, 4),
      PH_violated_OS_mv     = ph_violated_mv,
      HR_PFS                = round(hr_pfs,    4),
      CI95_low_PFS          = round(ci_low_pfs, 4),
      CI95_up_PFS           = round(ci_up_pfs,  4),
      P_PFS_LogRank         = p_logrank_pfs,
      P_PFS_Cox             = p_cox_pfs,
      PH_p_PFS              = round(ph_p_pfs, 4),
      PH_violated_PFS       = ph_violated_pfs,
      Limitations           = lim_str,
      stringsAsFactors      = FALSE
    )

    list(result = result_row, dist = dist_row, excl = excl_log_local)

  }, error = function(e) {
    # =========================================================
    # ERROR DIAGNOSIS
    # =========================================================

    err_call    <- if (!is.null(e$call)) deparse(e$call)[1] else "<unknown>"
    err_class   <- paste(class(e), collapse = ", ")
    err_time    <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

    ctx_lines <- if (length(step_context) > 0) {
      sapply(names(step_context), function(k) {
        val <- step_context[[k]]
        val_str <- tryCatch(
          paste(deparse(val, width.cutoff = 60), collapse = " "),
          error = function(e2) "<unprintable>"
        )
        sprintf("    %s = %s", k, substr(val_str, 1, 200))
      })
    } else "    <none>"

    message(sprintf(
      "  [%s] %s: ERROR at step '%s' [%s]: %s",
      project, gene, step_name, err_class, e$message
    ))

    err_file <- file.path(OUT_DIR, "diagnostics",
                          sprintf("%s_%s_ERROR.txt", project, gene))
    tryCatch({
      err_lines <- c(
        "================================================================",
        sprintf("ERROR REPORT for [%s / %s]", project, gene),
        "================================================================",
        sprintf("Time:           %s", err_time),
        sprintf("Failed step:    %s", step_name),
        sprintf("Error class:    %s", err_class),
        sprintf("Error message:  %s", e$message),
        sprintf("Failing call:   %s", err_call),
        "",
        "Step context (variable state at time of failure):",
        ctx_lines,
        "",
        "Session info:",
        sprintf("  R version:      %s", R.version.string),
        sprintf("  TCGAbiolinks:   %s", as.character(utils::packageVersion("TCGAbiolinks"))),
        sprintf("  survival:       %s", as.character(utils::packageVersion("survival"))),
        sprintf("  survminer:      %s", as.character(utils::packageVersion("survminer"))),
        sprintf("  dplyr:          %s", as.character(utils::packageVersion("dplyr"))),
        "================================================================"
      )
      writeLines(err_lines, err_file)
    }, error = function(e2) {
      message(sprintf("  [%s] %s: Failed to write error file: %s",
                      project, gene, e2$message))
    })

    tryCatch(while (dev.cur() > 1) dev.off(), error = function(e2) NULL)

    list(result = NULL, dist = NULL, excl = excl_log_local)
  })

  return(out)
}


# --- 6b. ------------------------------------------------------

analyze_project <- function(project) {

  message(sprintf("\n>>> [%s] Analyzing...", project))

  proj_step_name    <- "init"
  proj_step_context <- list()
  proj_set_step <- function(name, ctx = list()) {
    proj_step_name    <<- name
    proj_step_context <<- ctx
  }

  write_project_error <- function(e, stage_label) {
    err_call  <- if (!is.null(e$call)) deparse(e$call)[1] else "<unknown>"
    err_class <- paste(class(e), collapse = ", ")
    err_time  <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

    ctx_lines <- if (length(proj_step_context) > 0) {
      sapply(names(proj_step_context), function(k) {
        val_str <- tryCatch(
          paste(deparse(proj_step_context[[k]], width.cutoff = 60), collapse = " "),
          error = function(e2) "<unprintable>"
        )
        sprintf("    %s = %s", k, substr(val_str, 1, 200))
      })
    } else "    <none>"

    message(sprintf("[%s] %s ERROR at step '%s' [%s]: %s",
                    project, stage_label, proj_step_name, err_class, e$message))

    err_file <- file.path(OUT_DIR, "diagnostics",
                          sprintf("%s_PROJECT_%s_ERROR.txt", project, stage_label))
    tryCatch({
      writeLines(c(
        "================================================================",
        sprintf("PROJECT-LEVEL ERROR: %s [%s]", project, stage_label),
        "================================================================",
        sprintf("Time:           %s", err_time),
        sprintf("Failed step:    %s", proj_step_name),
        sprintf("Error class:    %s", err_class),
        sprintf("Error message:  %s", e$message),
        sprintf("Failing call:   %s", err_call),
        "",
        "Step context:",
        ctx_lines,
        "================================================================"
      ), err_file)
    }, error = function(e2) NULL)
  }

  loaded <- tryCatch({
    gdc_dir <- file.path(DOWNLOAD_BASE, project)

    proj_set_step("GDCquery", list(gdc_dir = gdc_dir))
    query <- GDCquery(project = project,
                      data.category = "Transcriptome Profiling",
                      data.type     = "Gene Expression Quantification",
                      workflow.type = "STAR - Counts")

    proj_set_step("GDCprepare")
    suppressMessages(se <- GDCprepare(query, directory = gdc_dir))

    proj_set_step("GDCquery_clinic")
    clin <- GDCquery_clinic(project = project, type = "clinical")

    proj_set_step("rowData_extraction",
                  list(n_clin = nrow(clin), n_se_rows = nrow(se),
                       assays = paste(assayNames(se), collapse = ",")))
    rd <- as.data.frame(rowData(se))

    list(se = se, clin = clin, rd = rd, n_raw_clin = nrow(clin))
  }, error = function(e) {
    write_project_error(e, "LOAD")
    NULL
  })

  if (is.null(loaded)) {
    return(list(results = NULL, distributions = NULL))
  }

  se         <- loaded$se
  clin       <- loaded$clin
  rd         <- loaded$rd
  n_raw_clin <- loaded$n_raw_clin

  tryCatch({
    diag_info <- data.frame(
      Project       = project,
      N_clinical    = n_raw_clin,
      N_genes       = nrow(rd),
      Has_age       = "age_at_diagnosis"      %in% colnames(clin),
      Has_stage     = "ajcc_pathologic_stage" %in% colnames(clin),
      Has_death     = "days_to_death"         %in% colnames(clin),
      Has_followup  = "days_to_last_follow_up" %in% colnames(clin),
      Has_vital     = "vital_status"          %in% colnames(clin),
      Has_recurrence= "days_to_recurrence"    %in% colnames(clin),
      Has_FPKM      = "fpkm_unstrand"         %in% assayNames(se),
      stringsAsFactors = FALSE
    )
    write_csv(diag_info,
              file.path(OUT_DIR, "diagnostics",
                        sprintf("%s_diagnostics.csv", project)))
  }, error = function(e) {
    message(sprintf("  [%s] diagnostics write failed: %s", project, e$message))
  })

  clin <- tryCatch({
    proj_set_step("clin_prep",
                  list(has_age   = "age_at_diagnosis"      %in% colnames(clin),
                       has_stage = "ajcc_pathologic_stage" %in% colnames(clin),
                       n_clin    = nrow(clin)))
    has_age   <- "age_at_diagnosis"      %in% colnames(clin)
    has_stage <- "ajcc_pathologic_stage" %in% colnames(clin)

    clin$age_at_diagnosis <- if (has_age) {
      suppressWarnings(as.numeric(clin$age_at_diagnosis))
    } else NA_real_

    clin$stage_num <- if (has_stage) {
      parse_stage(clin$ajcc_pathologic_stage)
    } else NA_real_

    clin
  }, error = function(e) {
    write_project_error(e, "CLIN_PREP")
    NULL
  })

  if (is.null(clin)) {
    return(list(results = NULL, distributions = NULL))
  }

  gene_outputs <- lapply(TARGET_GENES, function(gene) {
    analyze_gene(project, gene, se, clin, rd, n_raw_clin)
  })

  # Aggregation
  results_list <- lapply(gene_outputs, function(x) x$result)
  dist_list    <- lapply(gene_outputs, function(x) x$dist)
  excl_list    <- lapply(gene_outputs, function(x) x$excl)

  project_results       <- bind_rows(results_list)
  project_distributions <- bind_rows(dist_list)
  project_excl          <- bind_rows(excl_list)

  if (nrow(project_excl) > 0) {
    tryCatch({
      write_csv(project_excl,
                file.path(OUT_DIR, "exclusion_logs",
                          sprintf("%s_exclusion_log.csv", project)))
    }, error = function(e) {
      message(sprintf("  [%s] exclusion log write failed: %s", project, e$message))
    })
  }

  # Cleanup
  rm(se, clin, rd); gc()

  list(results = project_results, distributions = project_distributions)
}


# --- 7. Execute ---------------------------------------------------------------

all_projects <- grep("^TCGA-", getGDCprojects()$project_id, value = TRUE)
message(sprintf("\nFound %d TCGA projects.", length(all_projects)))

message("\n=== Phase 1: Sequential Downloads ===")
download_status     <- vapply(all_projects, download_project, logical(1))
projects_downloaded <- all_projects[download_status]
message(sprintf("\nDownloaded: %d/%d projects.",
                length(projects_downloaded), length(all_projects)))

if (length(projects_downloaded) == 0) {
  stop("FATAL: Keine Projekte heruntergeladen.")
}

message("\n=== Phase 2: Parallel Analysis ===")
n_cores   <- suppressWarnings(as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK")))
if (is.na(n_cores) || n_cores < 1) n_cores <- 4
n_workers <- min(n_cores, length(projects_downloaded))
plan(multicore, workers = n_workers)
message(sprintf("Using %d parallel workers.", n_workers))

raw_results <- future_lapply(projects_downloaded, analyze_project,
                              future.seed = TRUE)


# --- 8. Aggregation -----------------------------------------------------------

final_df <- bind_rows(lapply(raw_results, function(x)
                             if (!is.null(x)) x$results))
dist_df  <- bind_rows(lapply(raw_results, function(x)
                             if (!is.null(x)) x$distributions))

n_succ_projects <- sum(vapply(raw_results, function(x)
                              !is.null(x) && !is.null(x$results) && nrow(x$results) > 0,
                              logical(1)))

message(sprintf("\nProjects with at least one successful gene analysis: %d/%d",
                n_succ_projects, length(raw_results)))

if (nrow(final_df) == 0 || !"Gene" %in% colnames(final_df)) {
  message("\nFATAL: final_df ist leer.")
  message("Pruefe individuelle Logs in: ", file.path(OUT_DIR, "diagnostics"))
  stop("No successful entity-gene combinations.")
}

message(sprintf("Total entity-gene pairs: %d", nrow(final_df)))


final_df <- final_df %>%
  group_by(Gene) %>%
  mutate(
    FDR_OS_LogRank      = p.adjust(P_OS_LogRank,      method = "BH"),
    FDR_OS_univariate   = p.adjust(P_OS_univariate,   method = "BH"),
    FDR_OS_multivariate = p.adjust(P_OS_multivariate, method = "BH"),
    FDR_PFS_LogRank     = p.adjust(P_PFS_LogRank,     method = "BH"),
    Sig_OS_uv  = !is.na(FDR_OS_LogRank)      & FDR_OS_LogRank      < FDR_THRESHOLD,
    Sig_OS_mv  = !is.na(FDR_OS_multivariate) & FDR_OS_multivariate < FDR_THRESHOLD,
    Sig_PFS    = !is.na(FDR_PFS_LogRank)     & FDR_PFS_LogRank     < FDR_THRESHOLD,
    Direction_OS = dplyr::case_when(
      Sig_OS_uv & HR_OS_uv > 1 ~ "poor_prognosis",
      Sig_OS_uv & HR_OS_uv < 1 ~ "good_prognosis",
      TRUE                     ~ "ns"
    )
  ) %>%
  ungroup() %>%
  arrange(Gene, FDR_OS_LogRank)


# --- 9. Export ----------------------------------------------------------------

write_csv(final_df, file.path(OUT_DIR, "tables",
                               "MASTER_PanCancer_Survival_REMARK.csv"))
write_csv(dist_df,  file.path(OUT_DIR, "tables", "Marker_Distributions.csv"))

sig_df <- final_df %>% filter(Sig_OS_uv | Sig_OS_mv | Sig_PFS)
if (nrow(sig_df) > 0) {
  write_csv(sig_df, file.path(OUT_DIR, "tables", "Significant_Hits_FDR05.csv"))
}

ph_flags <- final_df %>%
  filter(isTRUE(PH_violated_OS_uv) |
         isTRUE(PH_violated_OS_mv) |
         isTRUE(PH_violated_PFS))
if (nrow(ph_flags) > 0) {
  write_csv(ph_flags, file.path(OUT_DIR, "tables", "PH_Violations.csv"))
}

excl_files <- list.files(file.path(OUT_DIR, "exclusion_logs"), full.names = TRUE)
if (length(excl_files) > 0) {
  all_excl <- bind_rows(lapply(excl_files, read_csv, show_col_types = FALSE))
  write_csv(all_excl,
            file.path(OUT_DIR, "tables", "Exclusion_Log_consolidated.csv"))
}


# --- 10. Optional Cleanup -----------------------------------------------------

if (CLEANUP_AFTER_RUN && dir.exists(DOWNLOAD_BASE)) {
  message(sprintf("\nCleaning up %s ...", DOWNLOAD_BASE))
  unlink(DOWNLOAD_BASE, recursive = TRUE)
}


# --- 11. Konsolidierte Fehlerübersicht ---------------------------------------

err_files <- list.files(file.path(OUT_DIR, "diagnostics"),
                         pattern = "_ERROR\\.txt$", full.names = TRUE)

if (length(err_files) > 0) {
  err_summary <- do.call(rbind, lapply(err_files, function(f) {
    lines <- tryCatch(readLines(f), error = function(e) character(0))
    extract <- function(prefix) {
      m <- grep(prefix, lines, value = TRUE)
      if (length(m) > 0) sub(paste0("^", prefix, "\\s*"), "", m[1]) else NA_character_
    }
    data.frame(
      file          = basename(f),
      project       = extract("Project|PROJECT-LEVEL ERROR:"),
      failed_step   = extract("Failed step:"),
      error_class   = extract("Error class:"),
      error_message = extract("Error message:"),
      stringsAsFactors = FALSE
    )
  }))
  write_csv(err_summary,
            file.path(OUT_DIR, "tables", "ERROR_SUMMARY.csv"))
  message(sprintf("\nERROR_SUMMARY.csv geschrieben (%d Fehler dokumentiert).",
                  nrow(err_summary)))
}


# --- 12. Final Report ---------------------------------------------------------

message(sprintf(
"\n================================================================
   REMARK Analysis complete
================================================================
   Total entity-gene pairs:        %d
   Significant OS  (FDR < %.2f):   %d (univariate)
   Significant OS  (FDR < %.2f):   %d (multivariate)
   Significant PFS (FDR < %.2f):   %d
   PH assumption violations:       %d
----------------------------------------------------------------
   Outputs:
     Master CSV:    %s/tables/MASTER_PanCancer_Survival_REMARK.csv
     KM plots:      %s/KM_plots/
     Raw plot CSVs: %s/raw_plot_data/
     Exclusions:    %s/tables/Exclusion_Log_consolidated.csv
     Diagnostics:   %s/diagnostics/
================================================================",
  nrow(final_df),
  FDR_THRESHOLD, sum(final_df$Sig_OS_uv,  na.rm = TRUE),
  FDR_THRESHOLD, sum(final_df$Sig_OS_mv,  na.rm = TRUE),
  FDR_THRESHOLD, sum(final_df$Sig_PFS,    na.rm = TRUE),
  nrow(ph_flags),
  OUT_DIR, OUT_DIR, OUT_DIR, OUT_DIR, OUT_DIR
))
