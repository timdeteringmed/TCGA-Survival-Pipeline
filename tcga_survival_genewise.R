#!/usr/bin/env Rscript

# ============================================
# \_\_\_\_  \_  \_\_\_\_\_  \_   Rheinisch-
# \_   \_\_ \__ \_  \_  \_  \_   Westfaelische
# \_\_\_  \_\_\_\_  \_  \_\_\_   Technische
# \_  \_   \__ \__  \_  \_  \_   Hochschule
# \_   \_   \_  \_  \_  \_  \_   Aachen
# ============================================

# ==============================================================================
# Author: Tim Detering
# Focus: TCGA genewise survival Pipeline
# filename: tcga_survival_genewise.R
# ==============================================================================

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(survival)
  library(survminer)
  library(future.apply)
  library(dplyr)
  library(readr)
})

# Parameter
TARGET_GENES <- c("ITIH2", "ITIH5")
OUT_DIR <- "results_tcga_final"
MIN_SAMPLES <- 30
dir.create(OUT_DIR, showWarnings = FALSE)

# HPC Setup
n_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
if (is.na(n_cores)) n_cores <- 4
plan(multicore, workers = n_cores)

analyze_project <- function(project) {
  project_results <- list()
  
  tryCatch({
    clin <- GDCquery_clinic(project = project, type = "clinical")
    query <- GDCquery(
      project = project,
      data.category = "Transcriptome Profiling",
      data.type = "Gene Expression Quantification",
      workflow.type = "STAR - Counts"
    )
    GDCdownload(query, method = "api", directory = "GDCdata")
    se <- GDCprepare(query, directory = "GDCdata")
    
    expr_matrix <- assay(se, "fpkm_unstrand")
    rownames(expr_matrix) <- rowData(se)$gene_name
    colnames(expr_matrix) <- substr(colnames(expr_matrix), 1, 12)
    expr_matrix <- t(simplify2array(by(t(expr_matrix), colnames(expr_matrix), colMeans)))

    # GENEWISE ITERATION
    for (gene in TARGET_GENES) {
      # Prüfung: Gen vorhanden?
      if (!(gene %in% rownames(expr_matrix))) {
        message(sprintf("[%s] SKIP %s: Gen nicht im Datensatz vorhanden.", project, gene))
        next
      }

      df_gene <- data.frame(
        submitter_id = colnames(expr_matrix),
        expression = expr_matrix[gene, ]
      )

      df_merged <- clin %>%
        inner_join(df_gene, by = "submitter_id") %>%
        mutate(
          OS_time = coalesce(as.numeric(days_to_death), as.numeric(days_to_last_follow_up)),
          OS_status = ifelse(vital_status == "Dead", 1, 0)
        ) %>%
        filter(!is.na(OS_time), OS_time > 0)

      if (nrow(df_merged) < MIN_SAMPLES) {
        message(sprintf("[%s] SKIP %s: Zu wenige Samples (N=%d) für KM-Plot.", project, gene, nrow(df_merged)))
        next
      }

      # OPTIMAL CUT OFF SELECTION WITH SURV_CUTPOINT
      res_cut <- try(surv_cutpoint(df_merged, time = "OS_time", event = "OS_status", variables = "expression"), silent = TRUE)
      if (inherits(res_cut, "try-error")) {
        message(sprintf("[%s] SKIP %s: Cut-Off konnte nicht berechnet werden.", project, gene))
        next
      }
      
      opt_cut <- res_cut$expression$estimate
      df_merged$Group <- ifelse(df_merged$expression >= opt_cut, "High", "Low")
      df_merged$Group <- factor(df_merged$Group, levels = c("Low", "High"))

      # KM-PLOTTING
      fit <- survfit(Surv(OS_time, OS_status) ~ Group, data = df_merged)
      p <- ggsurvplot(fit, data = df_merged, pval = TRUE, risk.table = TRUE,
                      title = paste(project, "-", gene), palette = c("blue", "red"))
      
      pdf(file.path(OUT_DIR, paste0(project, "_", gene, "_KM.pdf")), width = 7, height = 7, onefile = FALSE)
      print(p)
      dev.off()

      cox_fit <- coxph(Surv(OS_time, OS_status) ~ Group, data = df_merged)
      cox_sum <- summary(cox_fit)
      
      project_results[[gene]] <- data.frame(
        Entity = project,
        Gene = gene,
        N = nrow(df_merged),
        CutOff = opt_cut,
        HazardRatio = cox_sum$coefficients[1, "exp(coef)"],
        Lower_CI = cox_sum$conf.int[1, "lower .95"],
        Upper_CI = cox_sum$conf.int[1, "upper .95"],
        Cox_Pval = cox_sum$coefficients[1, "Pr(>|z|)"],
        KM_LogRank_Pval = surv_pvalue(fit)$pval
      )
    }
    
    rm(se, expr_matrix, query)
    gc()
    return(bind_rows(project_results))
    
  }, error = function(e) {
    message(sprintf("[%s] FATAL ERROR: %s", project, e$message))
    return(NULL)
  })
}

# Durchführung
tcga_projects <- grep("^TCGA-", getGDCprojects()$project_id, value = TRUE)
final_results <- future_lapply(tcga_projects, analyze_project, future.seed = TRUE)

# Export
all_stats <- bind_rows(final_results)
write_csv(all_stats, file.path(OUT_DIR, "TCGA_ITIH_Comprehensive_Results.csv"))
message("Analysis finished. Filen in: ", OUT_DIR)