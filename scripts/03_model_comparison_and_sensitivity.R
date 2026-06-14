# 03_model_comparison_and_sensitivity.R
#
# Supports:
# - Figure 2
# - Supplementary Table 1
# - Supplementary Table 3
# - Supplementary Figure 1
# - Response letter: GTEx lung FUSION TWAS sensitivity analysis

source("scripts/00_helpers.R")
paths <- read_paths()
out_dir <- make_outdir(file.path(paths[["out_dir"]], "03_model_comparison"))

compare_eas_eur_twas <- function(eas_file, eur_file, out_prefix) {
  eas <- read_table_auto(eas_file)
  eur <- read_table_auto(eur_file)

  if (!"z_select" %in% names(eas)) eas$z_select <- select_model_z(eas)
  if (!"z_select" %in% names(eur)) eur$z_select <- select_model_z(eur)
  if (!"fdr" %in% names(eas)) eas$fdr <- safe_fdr(eas$p_ACAT)
  if (!"fdr_ACAT" %in% names(eur)) eur$fdr_ACAT <- safe_fdr(eur$p_ACAT)

  eas2 <- eas %>%
    transmute(
      gene_id,
      gene_name = coalesce(!!!syms(intersect(c("gene_name", "gene_symbol"), names(eas))[1])),
      chr = chromosome,
      pos = gene_pos,
      z_EAS = z_select,
      p_EAS = p_ACAT,
      sig_EAS = fdr < 0.05
    )
  eur2 <- eur %>%
    transmute(
      gene_id,
      gene_name_EUR = coalesce(!!!syms(intersect(c("gene_name", "gene_symbol"), names(eur))[1])),
      z_EUR = z_select,
      p_EUR = p_ACAT,
      sig_EUR = fdr_ACAT < 0.05
    )

  both <- inner_join(eas2, eur2, by = "gene_id") %>%
    mutate(
      classification = case_when(
        sig_EAS & sig_EUR ~ "Shared",
        sig_EAS & !sig_EUR ~ "EAS-only",
        !sig_EAS & sig_EUR ~ "EUR-only",
        TRUE ~ "Not significant"
      )
    )

  significant_union <- full_join(eas2, eur2, by = "gene_id") %>%
    filter(coalesce(sig_EAS, FALSE) | coalesce(sig_EUR, FALSE)) %>%
    mutate(
      classification = case_when(
        coalesce(sig_EAS, FALSE) & coalesce(sig_EUR, FALSE) ~ "Shared",
        coalesce(sig_EAS, FALSE) & !coalesce(sig_EUR, FALSE) ~ "EAS-only",
        !coalesce(sig_EAS, FALSE) & coalesce(sig_EUR, FALSE) ~ "EUR-only",
        TRUE ~ "Not significant"
      )
    )

  write_tsv(significant_union, paste0(out_prefix, "_Supplementary_Table_3_EAS_vs_EUR_TWAS.tsv"))

  p <- ggplot() +
    geom_point(
      data = filter(both, classification == "Not significant"),
      aes(x = z_EUR, y = z_EAS),
      color = "grey75", alpha = 0.5, size = 0.8
    ) +
    geom_point(
      data = filter(both, classification != "Not significant"),
      aes(x = z_EUR, y = z_EAS, color = classification),
      alpha = 0.9, size = 1.5
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
    labs(x = "TWAS Z-score using EUR models", y = "TWAS Z-score using EAS models") +
    theme_pg()

  ggsave(paste0(out_prefix, "_Supplementary_Figure_1_EAS_vs_EUR_TWAS.pdf"), p, width = 6, height = 5)
  invisible(list(overlap = both, table_s3 = significant_union))
}

summarise_gtex_lung_twas <- function(gtex_file, out_prefix) {
  gtex <- read_table_auto(gtex_file)
  p_col <- intersect(c("P", "PANEL.P", "pvalue", "p_ACAT"), names(gtex))[1]
  if (is.na(p_col)) stop("GTEx lung TWAS file must include a P-value column.")
  gtex$fdr <- safe_fdr(gtex[[p_col]])
  write_tsv(gtex, paste0(out_prefix, "_GTEx_lung_TWAS_processed.tsv"))
  invisible(gtex)
}

out_prefix <- file.path(out_dir, "model_comparison")
compare_eas_eur_twas(paths[["eas_twas_njmu"]], paths[["eur_twas_njmu"]], out_prefix)
summarise_gtex_lung_twas(paths[["gtex_lung_twas_njmu"]], out_prefix)
