# Shared helper functions for the EAS lung cancer TWAS/MWAS workflow.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

read_paths <- function(path_file = "config/paths_template.tsv") {
  paths <- fread(path_file)
  stopifnot(all(c("key", "path") %in% names(paths)))
  x <- setNames(paths$path, paths$key)
  return(x)
}

make_outdir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  invisible(path)
}

read_table_auto <- function(path) {
  if (grepl("\\.csv(\\.gz)?$", path, ignore.case = TRUE)) {
    fread(path)
  } else {
    fread(path, sep = "\t")
  }
}

write_tsv <- function(x, path) {
  make_outdir(dirname(path))
  fwrite(x, path, sep = "\t", na = "NA")
}

safe_fdr <- function(p) {
  out <- rep(NA_real_, length(p))
  keep <- !is.na(p)
  out[keep] <- p.adjust(p[keep], method = "fdr")
  out
}

select_model_z <- function(df, model_col = "model_best") {
  models <- as.character(df[[model_col]])
  z <- rep(NA_real_, nrow(df))
  for (i in seq_len(nrow(df))) {
    z_col <- paste0("z_", models[i])
    if (z_col %in% names(df)) z[i] <- suppressWarnings(as.numeric(df[[z_col]][i]))
  }
  z
}

acat <- function(p, weights = NULL) {
  p <- as.numeric(p)
  p <- p[!is.na(p)]
  if (!length(p)) return(NA_real_)
  if (any(p < 0 | p > 1)) stop("ACAT received P values outside [0, 1].")
  if (is.null(weights)) weights <- rep(1 / length(p), length(p))
  weights <- weights / sum(weights)
  if (any(p == 0)) return(0)
  if (any(p == 1)) p[p == 1] <- 1 - .Machine$double.eps
  stat <- sum(weights * tan((0.5 - p) * pi))
  pcauchy(stat, lower.tail = FALSE)
}

classify_discovery_replication <- function(discovery_p, replication_p,
                                           discovery_fdr = NULL,
                                           fdr_cutoff = 0.05,
                                           replication_p_cutoff = 0.05) {
  discovery_sig <- if (is.null(discovery_fdr)) {
    safe_fdr(discovery_p) < fdr_cutoff
  } else {
    discovery_fdr < fdr_cutoff
  }
  replicated <- discovery_sig & replication_p < replication_p_cutoff
  list(discovery_sig = discovery_sig, replicated = replicated)
}

theme_pg <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text = element_text(color = "black"),
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0)
    )
}
