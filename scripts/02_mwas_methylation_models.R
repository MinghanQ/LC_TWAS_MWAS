# 02_mwas_methylation_models.R
#
# Supports:
# - Methods: DNA methylation prediction model construction/application
# - Table 2
# - Supplementary Table 4
# - Figure 3C-D
#
# The methylation model weights and covariance files are shared as data. This
# script contains the cleaned MetaMeth scan logic and downstream MWAS summaries.

source("scripts/00_helpers.R")
paths <- read_paths()
out_dir <- make_outdir(file.path(paths[["out_dir"]], "02_mwas"))

metameth_one <- function(gwas_data, model_data, ld_mat, n = NULL,
                         snp_id = "SNP", lambda = 0.1) {
  if (!snp_id %in% names(model_data)) stop("Model data lacks SNP column: ", snp_id)
  model_data <- as.data.table(model_data)
  gwas_data <- as.data.table(gwas_data)
  setkeyv(model_data, snp_id)
  setkeyv(gwas_data, snp_id)
  dat <- merge(model_data, gwas_data)

  result <- list(nsnps_model = nrow(model_data), nsnps_used = nrow(dat),
                 z = NA_real_, p_norm = NA_real_, p_student = NA_real_)
  if (!nrow(dat)) return(result)

  keep <- dat$GWAS_A2 == dat$Allele2 | dat$GWAS_A2 == dat$Allele0
  dat <- dat[keep]
  if (!nrow(dat)) return(result)

  flip <- dat$GWAS_A2 != dat$Allele2
  dat$GWAS_Beta[flip] <- -dat$GWAS_Beta[flip]

  colnames(ld_mat) <- gsub("\\.", ":", gsub("snp_", "", colnames(ld_mat)))
  idx <- match(dat[[snp_id]], colnames(ld_mat))
  dat <- dat[!is.na(idx)]
  idx <- idx[!is.na(idx)]
  if (!nrow(dat)) return(result)

  gamma <- as.matrix(ld_mat[idx, idx])
  sigma_l <- sqrt(diag(gamma))
  diag(gamma) <- diag(gamma) + lambda * diag(gamma)

  weights <- matrix(dat$Beta, ncol = 1)
  numerator <- sum(weights * sigma_l * dat$GWAS_Beta)
  denominator <- as.numeric(t(weights) %*% gamma %*% weights)
  z <- numerator / sqrt(denominator)

  result$nsnps_used <- nrow(dat)
  result$z <- z
  result$p_norm <- 2 * pnorm(abs(z), lower.tail = FALSE)
  if (!is.null(n)) result$p_student <- 2 * pt(abs(z), df = n - 2, lower.tail = FALSE)
  result
}

process_mwas_results <- function(njmu_file, bbj_file, meta_file, out_prefix) {
  njmu <- read_table_auto(njmu_file)
  bbj <- read_table_auto(bbj_file)
  meta <- read_table_auto(meta_file)

  p_col <- intersect(c("Pnorm", "P", "P.value", "Pmeta", "p_meta"), names(njmu))[1]
  if (!is.na(p_col) && !"fdr" %in% names(njmu)) {
    njmu$fdr <- safe_fdr(njmu[[p_col]])
  }

  write_tsv(njmu, paste0(out_prefix, "_MWAS_NJMU_processed.tsv"))
  write_tsv(bbj, paste0(out_prefix, "_MWAS_BBJ_processed.tsv"))
  write_tsv(meta, paste0(out_prefix, "_Supplementary_Table_4_MWAS.tsv"))

  if (all(c("CHR_hg38", "Pos_C_hg38") %in% names(meta))) {
    p_meta_col <- intersect(c("Pmeta", "P.value", "Pnorm", "p_meta"), names(meta))[1]
    mh <- meta %>%
      mutate(
        chr = as.factor(CHR_hg38),
        pos = as.numeric(Pos_C_hg38),
        logp = pmin(-log10(pmax(as.numeric(.data[[p_meta_col]]), .Machine$double.xmin)), 50)
      )
    p_mh <- ggplot(mh, aes(x = pos, y = logp)) +
      geom_point(size = 0.6, alpha = 0.75) +
      facet_grid(. ~ chr, scales = "free_x", space = "free_x") +
      labs(x = "Genomic position", y = expression(-log[10](P))) +
      theme_pg(8) +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    ggsave(paste0(out_prefix, "_Figure3C_MWAS_manhattan.pdf"), p_mh, width = 9, height = 3.5)
  }

  invisible(list(njmu = njmu, bbj = bbj, meta = meta))
}

out_prefix <- file.path(out_dir, "mwas")
process_mwas_results(
  paths[["mwas_njmu"]],
  paths[["mwas_bbj"]],
  paths[["mwas_meta"]],
  out_prefix
)
