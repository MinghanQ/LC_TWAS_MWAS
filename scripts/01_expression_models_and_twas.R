# 01_expression_models_and_twas.R
#
# Supports:
# - Methods: EAS expression prediction model construction/evaluation
# - Figure 2
# - Table 1
# - Supplementary Table 1
# - Supplementary Table 2
#
# Notes:
# - The full SUMMIT model training requires external SUMMIT source files and
#   large eQTL/model inputs. Public model weights should be shared separately.
# - This release script standardizes model-performance and TWAS summary outputs
#   used in the manuscript.

source("scripts/00_helpers.R")
paths <- read_paths()
out_dir <- make_outdir(file.path(paths[["out_dir"]], "01_expression_twas"))

summarise_expression_model_performance <- function(eas_perf, eur_perf, out_prefix) {
  eas <- read_table_auto(eas_perf) %>% mutate(ancestry = "EAS")
  eur <- read_table_auto(eur_perf) %>% mutate(ancestry = "EUR")
  perf <- bind_rows(eas, eur)

  r2_col <- intersect(c("r2_test", "R2", "test_R2", "r2"), names(perf))[1]
  method_col <- intersect(c("model", "method", "model_best", "Method"), names(perf))[1]
  if (is.na(r2_col) || is.na(method_col)) {
    stop("Model performance files must include an R2 column and a model/method column.")
  }

  perf <- perf %>%
    mutate(
      r2_value = as.numeric(.data[[r2_col]]),
      method = as.character(.data[[method_col]]),
      qualified = r2_value > 0.01
    )

  table_s1 <- perf %>%
    group_by(method, ancestry) %>%
    summarise(
      n_genes_r2_gt_0_01 = sum(qualified, na.rm = TRUE),
      median_r2 = median(r2_value[qualified], na.rm = TRUE),
      .groups = "drop"
    )

  write_tsv(table_s1, paste0(out_prefix, "_Supplementary_Table_1_model_performance.tsv"))

  p_count <- ggplot(table_s1, aes(x = method, y = n_genes_r2_gt_0_01, fill = ancestry)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.65) +
    labs(x = NULL, y = "Genes with test R2 > 0.01") +
    theme_pg()

  p_r2 <- ggplot(filter(perf, qualified), aes(x = ancestry, y = r2_value, fill = ancestry)) +
    geom_boxplot(width = 0.55, outlier.size = 0.35) +
    facet_wrap(~ method, scales = "free_y") +
    labs(x = NULL, y = "Test R2") +
    theme_pg()

  ggsave(paste0(out_prefix, "_Figure2A_model_count.pdf"), p_count, width = 6.2, height = 4.2)
  ggsave(paste0(out_prefix, "_Figure2B_model_r2.pdf"), p_r2, width = 7.5, height = 4.5)

  invisible(list(performance = perf, table_s1 = table_s1))
}

process_twas_results <- function(njmu_file, bbj_file, meta_file, out_prefix) {
  njmu <- read_table_auto(njmu_file)
  bbj <- read_table_auto(bbj_file)
  meta <- read_table_auto(meta_file)

  if (!"z_best" %in% names(njmu) && "model_best" %in% names(njmu)) {
    njmu$z_best <- select_model_z(njmu)
  }
  if (!"z_best" %in% names(bbj) && "model_best" %in% names(bbj)) {
    bbj$z_best <- select_model_z(bbj)
  }

  if (!"fdr" %in% names(njmu) && "p_ACAT" %in% names(njmu)) {
    njmu$fdr <- safe_fdr(njmu$p_ACAT)
  }

  table_s2 <- meta %>%
    mutate(
      p_meta = coalesce(!!!syms(intersect(c("p_select", "Pmeta", "p_meta", "p_ACAT"), names(meta))[1])),
      z_meta = coalesce(!!!syms(intersect(c("z_select", "Zmeta", "z_meta"), names(meta))[1]))
    )

  write_tsv(njmu, paste0(out_prefix, "_TWAS_NJMU_processed.tsv"))
  write_tsv(bbj, paste0(out_prefix, "_TWAS_BBJ_processed.tsv"))
  write_tsv(table_s2, paste0(out_prefix, "_Supplementary_Table_2_TWAS.tsv"))

  if (all(c("chromosome", "gene_pos", "p_ACAT") %in% names(njmu))) {
    mh <- njmu %>%
      mutate(
        chromosome = as.factor(chromosome),
        logp = pmin(-log10(pmax(p_ACAT, .Machine$double.xmin)), 50)
      )
    p_mh <- ggplot(mh, aes(x = gene_pos, y = logp)) +
      geom_point(size = 0.6, alpha = 0.75) +
      facet_grid(. ~ chromosome, scales = "free_x", space = "free_x") +
      labs(x = "Genomic position", y = expression(-log[10](P[ACAT]))) +
      theme_pg(8) +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    ggsave(paste0(out_prefix, "_Figure3A_TWAS_manhattan.pdf"), p_mh, width = 9, height = 3.5)
  }

  invisible(list(njmu = njmu, bbj = bbj, meta = table_s2))
}

out_prefix <- file.path(out_dir, "expression_twas")
summarise_expression_model_performance(
  paths[["eas_expression_model_performance"]],
  paths[["eur_expression_model_performance"]],
  out_prefix
)

process_twas_results(
  paths[["eas_twas_njmu"]],
  paths[["eas_twas_bbj"]],
  paths[["eas_twas_meta"]],
  out_prefix
)
