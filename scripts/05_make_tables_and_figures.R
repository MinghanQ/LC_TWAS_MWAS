# 05_make_tables_and_figures.R
#
# Final assembly script for manuscript-facing outputs.
#
# Supports:
# - Table 1
# - Table 2
# - Supplementary Tables 1-7
# - Figure 2
# - Figure 3
# - Figure 4 and Supplementary Figures 1-4 through processed inputs

source("scripts/00_helpers.R")
paths <- read_paths()
out_dir <- make_outdir(file.path(paths[["out_dir"]], "05_final_outputs"))

format_p <- function(p, digits = 3) {
  ifelse(is.na(p), "NA", formatC(as.numeric(p), format = "e", digits = digits))
}

make_validated_twas_table <- function(twas_meta, coloc_eqtl, smr_twas, out_file) {
  twas <- read_table_auto(twas_meta)
  coloc <- read_table_auto(coloc_eqtl)
  smr <- read_table_auto(smr_twas)

  gene_col <- intersect(c("gene_id", "MarkerName", "Gene"), names(twas))[1]
  table1 <- twas %>%
    left_join(coloc, by = setNames("MarkerName", gene_col)) %>%
    left_join(smr, by = setNames("Gene", gene_col))

  write_tsv(table1, out_file)
  invisible(table1)
}

make_validated_mwas_table <- function(mwas_meta, coloc_mqtl, smr_mwas, out_file) {
  mwas <- read_table_auto(mwas_meta)
  coloc <- read_table_auto(coloc_mqtl)
  smr <- read_table_auto(smr_mwas)

  cpg_col <- intersect(c("CpG", "gid", "MarkerName"), names(mwas))[1]
  table2 <- mwas %>%
    left_join(coloc, by = setNames("CpG", cpg_col)) %>%
    left_join(smr, by = setNames("CpG", cpg_col))

  write_tsv(table2, out_file)
  invisible(table2)
}

make_twas_mwas_z_scatter <- function(twas_njmu, twas_bbj, mwas_njmu, mwas_bbj, out_prefix) {
  twas1 <- read_table_auto(twas_njmu)
  twas2 <- read_table_auto(twas_bbj)
  mwas1 <- read_table_auto(mwas_njmu)
  mwas2 <- read_table_auto(mwas_bbj)

  if (!"z_best" %in% names(twas1)) twas1$z_best <- select_model_z(twas1)
  if (!"z_best" %in% names(twas2)) twas2$z_best <- select_model_z(twas2)

  twas <- inner_join(
    twas1 %>% transmute(gene_id, z_discovery = z_best),
    twas2 %>% transmute(gene_id, z_replication = z_best),
    by = "gene_id"
  )

  p_twas <- ggplot(twas, aes(x = z_discovery, y = z_replication)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
    geom_point(size = 1.2, alpha = 0.75) +
    labs(x = "Discovery TWAS Z-score", y = "Replication TWAS Z-score") +
    theme_pg()
  ggsave(paste0(out_prefix, "_Figure3B_TWAS_Z_scatter.pdf"), p_twas, width = 5.2, height = 4.5)

  z_col1 <- intersect(c("Zscore", "Z", "Zscore.njmu"), names(mwas1))[1]
  z_col2 <- intersect(c("Zscore", "Z", "Zscore.bbj"), names(mwas2))[1]
  cpg_col1 <- intersect(c("gid", "CpG", "MarkerName"), names(mwas1))[1]
  cpg_col2 <- intersect(c("gid", "CpG", "MarkerName"), names(mwas2))[1]

  mwas <- inner_join(
    mwas1 %>% transmute(CpG = .data[[cpg_col1]], z_discovery = as.numeric(.data[[z_col1]])),
    mwas2 %>% transmute(CpG = .data[[cpg_col2]], z_replication = as.numeric(.data[[z_col2]])),
    by = "CpG"
  )

  p_mwas <- ggplot(mwas, aes(x = z_discovery, y = z_replication)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
    geom_point(size = 1.2, alpha = 0.75) +
    labs(x = "Discovery MWAS Z-score", y = "Replication MWAS Z-score") +
    theme_pg()
  ggsave(paste0(out_prefix, "_Figure3D_MWAS_Z_scatter.pdf"), p_mwas, width = 5.2, height = 4.5)
}

out_prefix <- file.path(out_dir, "final")

make_validated_twas_table(
  paths[["eas_twas_meta"]],
  paths[["coloc_eqtl_gwas"]],
  paths[["smr_twas"]],
  paste0(out_prefix, "_Table1_validated_TWAS_genes.tsv")
)

make_validated_mwas_table(
  paths[["mwas_meta"]],
  paths[["coloc_mqtl_gwas"]],
  paths[["smr_mwas"]],
  paste0(out_prefix, "_Table2_validated_MWAS_CpGs.tsv")
)

make_twas_mwas_z_scatter(
  paths[["eas_twas_njmu"]],
  paths[["eas_twas_bbj"]],
  paths[["mwas_njmu"]],
  paths[["mwas_bbj"]],
  out_prefix
)
