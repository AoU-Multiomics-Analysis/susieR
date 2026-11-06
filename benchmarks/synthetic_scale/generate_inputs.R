suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 1L)
output_dir <- args[[1]]
n_samples <- if (length(args) >= 2L) as.integer(args[[2]]) else 10000L
n_variants <- if (length(args) >= 3L) as.integer(args[[3]]) else 20000L
stopifnot(n_samples >= 100L, n_variants >= 100L)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(20261105)
sample_ids <- as.character(100000L + seq_len(n_samples))
positions <- 100000L + seq_len(n_variants) * 50L
variant_ids <- paste0("chr1_", positions, "_A_G")
block_size <- 100L
n_blocks <- ceiling(n_variants / block_size)
message("Generate correlated diploid genotypes: ", n_samples, " samples, ", n_variants, " variants")
genotypes <- matrix(0L, nrow = n_variants, ncol = n_samples)
for (block in seq_len(n_blocks)) {
  rows <- seq.int((block - 1L) * block_size + 1L, min(block * block_size, n_variants))
  frequency <- runif(1L, 0.1, 0.45)
  haplotype_1 <- rbinom(n_samples, 1L, frequency)
  haplotype_2 <- rbinom(n_samples, 1L, frequency)
  for (row in rows) {
    flip_1 <- rbinom(n_samples, 1L, 0.01)
    flip_2 <- rbinom(n_samples, 1L, 0.01)
    genotypes[row, ] <- abs(haplotype_1 - flip_1) + abs(haplotype_2 - flip_2)
  }
}
covariates <- matrix(rnorm(n_samples * 20L), nrow = n_samples)
causal_rows <- unique(as.integer(round(c(0.15, 0.5, 0.85) * n_variants)))
effects <- c(0.4, -0.35, 0.3)[seq_along(causal_rows)]
phenotype <- as.numeric(crossprod(effects, genotypes[causal_rows, , drop = FALSE])) +
  0.2 * covariates[, 1L] - 0.15 * covariates[, 2L] + rnorm(n_samples)
missing_indices <- sample.int(length(genotypes), as.integer(length(genotypes) * 0.001))
midpoint <- length(missing_indices) %/% 2L
genotypes[missing_indices[seq_len(midpoint)]] <- -1L
genotypes[missing_indices[seq.int(midpoint + 1L, length(missing_indices))]] <- NA_integer_
colnames(genotypes) <- sample_ids
dosage_tsv <- file.path(output_dir, "dosages.tsv")
message("Write dosage table")
dosage <- data.table::as.data.table(genotypes)
rm(genotypes)
invisible(gc())
dosage <- data.table::data.table(CHROM = "1", POS = positions, REF = "A", ALT = "G") |>
  cbind(dosage)
data.table::fwrite(dosage, dosage_tsv, sep = "\t", quote = FALSE, na = "NA", nThread = 1L)
rm(dosage)
invisible(gc())
message("Compress and index dosage table")
dosage_bgz <- Rsamtools::bgzip(dosage_tsv, dest = paste0(dosage_tsv, ".bgz"), overwrite = TRUE)
Rsamtools::indexTabix(dosage_bgz, seq = 1L, start = 2L, end = 2L, skip = 1L, zeroBased = FALSE)
unlink(dosage_tsv)
phenotype_table <- tibble(`#chr` = "1", start = 999999L, end = 1000000L,
                          phenotype_id = c("synthetic_gene", "skip")) |>
  bind_cols(as_tibble(matrix(rep(phenotype, each = 2L), nrow = 2L,
                            dimnames = list(NULL, sample_ids))))
write_tsv(phenotype_table, file.path(output_dir, "phenotypes.bed"))
covariate_table <- as_tibble(t(covariates), .name_repair = "minimal")
names(covariate_table) <- sample_ids
covariate_table <- tibble(ID = paste0("covariate_", seq_len(ncol(covariates)))) |>
  bind_cols(covariate_table)
write_tsv(covariate_table, file.path(output_dir, "covariates.tsv"))
write_tsv(tibble(sample_id = sample_ids), file.path(output_dir, "samples.tsv"))
write_tsv(tibble(phenotype_id = "synthetic_gene", qval = 0.001), file.path(output_dir, "permutations.tsv"))
write_tsv(tibble(variant_id = variant_ids[causal_rows], effect = effects), file.path(output_dir, "causal_variants.tsv"))
jsonlite::write_json(list(samples = n_samples, variants = n_variants, covariates = 20L,
                          ld_block_size = block_size, haplotype_flip_probability = 0.01,
                          missing_fraction = 0.001, seed = 20261105L,
                          note = "One phenotype, three causal effects, correlated diploid genotypes. Half of missing dosages use -1; half use NA."),
                     file.path(output_dir, "fixture.json"), auto_unbox = TRUE, pretty = TRUE)
message("Synthetic inputs are ready")
