set.seed(20260827)

fixture_dir <- "tests/fixtures/checkpointed_window"
dir.create(fixture_dir, recursive = TRUE, showWarnings = FALSE)

sample_ids <- sprintf("sample_%02d", seq_len(40L))
genotype <- rbind(
  rep(c(0, 1, 2, 0, 1), 8L),
  rep(c(0, 0, 1, 2, 1, 0, 2, 1), 5L),
  rep(c(2, 1, 0, 1, 0), 8L),
  rep(c(0, 1, 0, 2), 10L),
  rep(c(2, 0, 1, 0, 0, 1, 2, 1), 5L),
  rep(c(0, 0, 0, 1, 1, 2, 0, 1, 2, 1), 4L)
)
genotype[3L, 4L] <- -1

dosage <- tibble::tibble(
  `#CHROM` = rep("1", 6L),
  POS = seq(100000L, 600000L, by = 100000L),
  REF = c("A", "C", "G", "T", "A", "C"),
  ALT = c("G", "T", "A", "C", "C", "G")
) |>
  dplyr::bind_cols(tibble::as_tibble(genotype, .name_repair = ~ sample_ids))
readr::write_tsv(dosage, file.path(fixture_dir, "window_dosage.tsv"))

shared_pc1 <- seq(-1, 1, length.out = length(sample_ids))
expression_pc1 <- sin(seq_along(sample_ids) / 4)
genotype_for_traits <- genotype
genotype_for_traits[genotype_for_traits == -1] <- NA_real_
genotype_for_traits[3L, 4L] <- mean(genotype_for_traits[3L, ], na.rm = TRUE)

expression_trait <-
  4.5 * genotype_for_traits[2L, ] -
  4.0 * genotype_for_traits[5L, ] +
  0.8 * shared_pc1 +
  0.6 * expression_pc1 +
  stats::rnorm(length(sample_ids), sd = 0.2)
splicing_trait <-
  3.5 * genotype_for_traits[4L, ] +
  0.5 * shared_pc1 +
  stats::rnorm(length(sample_ids), sd = 0.3)

phenotypes <- tibble::tibble(
  chromosome = c("12", "7"),
  start = c(1200000L, 700000L),
  end = c(1200001L, 700001L),
  phenotype_id = c("linked_splicing", "linked_expression")
) |>
  dplyr::bind_cols(
    tibble::as_tibble(
      rbind(splicing_trait, expression_trait),
      .name_repair = ~ sample_ids
    )
  )
readr::write_tsv(
  phenotypes,
  file.path(fixture_dir, "window_phenotypes.bed.gz")
)

manifest <- tibble::tribble(
  ~window_id, ~phenotype_id, ~modality, ~phenotype_file, ~p_value,
  "chr1_0_2000000", "linked_splicing", "splicing", "window_phenotypes.bed.gz", 2e-8,
  "chr1_0_2000000", "linked_expression", "expression", "window_phenotypes.bed.gz", 1e-12
)
readr::write_tsv(manifest, file.path(fixture_dir, "window_phenotypes.tsv"))

shared_covariates <- tibble::tibble(covariate_id = "shared_pc1") |>
  dplyr::bind_cols(tibble::as_tibble(t(shared_pc1), .name_repair = ~ sample_ids))
readr::write_tsv(
  shared_covariates,
  file.path(fixture_dir, "shared_covariates.tsv")
)

expression_covariates <- tibble::tibble(covariate_id = "expression_pc1") |>
  dplyr::bind_cols(tibble::as_tibble(t(expression_pc1), .name_repair = ~ sample_ids))
readr::write_tsv(
  expression_covariates,
  file.path(fixture_dir, "expression_covariates.tsv")
)

keep_samples <- tibble::tibble(sample_id = rev(sample_ids[5:40]))
readr::write_tsv(keep_samples, file.path(fixture_dir, "keep_samples.tsv"))

message("Created deterministic checkpointed-window fixtures.")
