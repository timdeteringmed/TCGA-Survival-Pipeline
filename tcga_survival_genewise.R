# ============================================ #
#  \_\_\_\_  \_  \_\_\_\_\_  \_  Rheinisch-    #
#  \_   \_\_ \__ \_  \_  \_  \_  Westfaelische #
#  \_\_\_  \_\_\_\_  \_  \_\_\_  Technische    #
#  \_  \_   \__ \__  \_  \_  \_  Hochschule    #
#  \_   \_   \_  \_  \_  \_  \_  Aachen        #
# ============================================ #

# ==============================================================================
# TCGA-Survival Analysis
# Project: MD
# Author: Tim Detering
# filename: TCGA_surv_analysis.R
# ==============================================================================

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(survival)
  library(survminer)
  library(future.apply)
  library(dplyr)
  library(readr)
})

# --- Parameter ---
TARGET_GENES <- c("ITIH2", "ITIH5")
OUT_DIR <- "results_tcga_final"
MIN_SAMPLES <- 30
dir.create(OUT_DIR, showWarnings = FALSE)

# HPC Setup (RWTH Cluster)
n_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
if (is.na(n_cores)) n_cores <- 4
plan(multicore, workers = n_cores)

analyze_project <- function(project) {
  project_results <- list()
  message(sprintf("\n[%s] Processing entity...", project))
  
  tryCatch({
    # 1. Daten-Download & Vorbereitung
    query <- GDCquery(
      project = project,
      data.category = "Transcriptome Profiling",
      data.type = "Gene Expression Quantification",
      workflow.type = "STAR - Counts"
    )
    GDCdownload(query, method = "api", directory = "GDCdata")
    se <- GDCprepare(query, directory = "GDCdata")
    
    clin <- GDCquery_clinic(project = project, type = "clinical")
    
    rd <- as.data.frame(rowData(se))
    
    for (gene in TARGET_GENES) {
      
      gene_idx <- which(rd$gene_name == gene)
      
      if (length(gene_idx) == 0) {
        message(sprintf("[%s] SKIP %s: Gene symbol not found in rowData.", project, gene))
        next
      }
      
      gene_idx <- gene_idx[1]
      
      expr_vals <- assay(se, "fpkm_unstrand")[gene_idx, ]
      
      patients <- substr(names(expr_vals), 1, 12)
      df_gene <- data.frame(submitter_id = patients, expression = expr_vals) %>%
        group_by(submitter_id) %>%
        summarise(expression = mean(expression), .groups = 'drop')
      
      df_merged <- clin %>%
        inner_join(df_gene, by = "submitter_id") %>%
        mutate(
          OS_time = coalesce(as.numeric(days_to_death), as.numeric(days_to_last_follow_up)),
          OS_status = ifelse(vital_status == "Dead", 1, 0)
        ) %>%
        filter(!is.na(OS_time), OS_time > 0)
      
      if (nrow(df_merged) < MIN_SAMPLES) {
        message(sprintf("[%s] SKIP %s: Too few samples with clinical data (N=%d).", project, gene, nrow(df_merged)))
        next
      }
      
      res_cut <- try(surv_cutpoint(df_merged, time = "OS_time", event = "OS_status", variables = "expression"), silent = TRUE)
      
      if (inherits(res_cut, "try-error")) {
        message(sprintf("[%s] SKIP %s: Could not determine optimal cut-point (check variance).", project, gene))
        next
      }
      
      opt_cut <- res_cut$expression$estimate
      df_merged$Group <- ifelse(df_merged$expression >= opt_cut, "High", "Low")
      df_merged$Group <- factor(df_merged$Group, levels = c("Low", "High"))
      
      fit <- survfit(Surv(OS_time, OS_status) ~ Group, data = df_merged)
      
      p <- ggsurvplot(fit, data = df_merged, pval = TRUE, risk.table = TRUE,
                      title = paste(project, "-", gene),
                      palette = c("blue", "red"), 
                      legend.labs = c("Low Expression", "High Expression"))
      
      pdf(file.path(OUT_DIR, paste0(project, "_", gene, "_KM.pdf")), width = 8, height = 7, onefile = FALSE)
      print(p)
      dev.off()
      
      cox_fit <- coxph(Surv(OS_time, OS_status) ~ Group, data = df_merged)
      cox_sum <- summary(cox_fit)
      
      project_results[[gene]] <- data.frame(
        Entity = project,
        Gene = gene,
        N = nrow(df_merged),
        CutOff_FPKM = opt_cut,
        HazardRatio = cox_sum$coefficients[1, "exp(coef)"],
        Cox_Pval = cox_sum$coefficients[1, "Pr(>|z|)"],
        KM_Pval = surv_pvalue(fit)$pval
      )
    }
    
    rm(se, query, clin, rd)
    gc()
    
    return(bind_rows(project_results))
    
  }, error = function(e) {
    message(sprintf("[%s] CRITICAL ERROR: %s", project, e$message))
    return(NULL)
  })
}

all_projects <- grep("^TCGA-", getGDCprojects()$project_id, value = TRUE)
master_results_list <- future_lapply(all_projects, analyze_project, future.seed = TRUE)

all_stats <- bind_rows(master_results_list)
write_csv(all_stats, file.path(OUT_DIR, "TCGA_ITIH_PanCancer_Summary.csv"))

message("\nDone! All results are in: ", OUT_DIR)
