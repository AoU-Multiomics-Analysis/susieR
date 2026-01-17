InitEmptyVariantDf <- function() {
#Define empty data frames
empty_variant_df = dplyr::tibble(
  molecular_trait_id = character(),
  variant = character(),
  chromosome = character(),
  position = integer(),
  ref = character(),
  alt = character(),
  cs_id = character(),
  cs_index = character(),
  low_purity = character(),
  region = character(),
  pip = numeric(),
  z = numeric(),
  posterior_mean = numeric(),
  posterior_sd = numeric(),
  X_column_scale_factors = numeric(),
  alpha1 = numeric(),
  alpha2 = numeric(),
  alpha3 = numeric(),
  alpha4 = numeric(),
  alpha5 = numeric(),
  alpha6 = numeric(),
  alpha7 = numeric(),
  alpha8 = numeric(),
  alpha9 = numeric(),
  alpha10 = numeric(),
  mu_1 = numeric(),
  mu_2 = numeric(),
  mu_3 = numeric(),
  mu_4 = numeric(),
  mu_5 = numeric(),
  mu_6 = numeric(),
  mu_7 = numeric(),
  mu_8 = numeric(),
  mu_9 = numeric(),
  mu_10 = numeric(),
  mu2_1 = numeric(),
  mu2_2 = numeric(),
  mu2_3 = numeric(),
  mu2_4 = numeric(),
  mu2_5 = numeric(),
  mu2_6 = numeric(),
  mu2_7 = numeric(),
  mu2_8 = numeric(),
  mu2_9 = numeric(),
  mu2_10 = numeric()
)
empty_variant_df


}


InitEmptyLbfDf <- function() {
empty_lbf_df = dplyr::tibble(
  molecular_trait_id = character(),
  region = character(),
  variant = character(),
  chromosome = character(),
  position = integer(),
  lbf_variable1 = numeric(),
  lbf_variable2 = numeric(),
  lbf_variable3 = numeric(),
  lbf_variable4 = numeric(),
  lbf_variable5 = numeric(),
  lbf_variable6 = numeric(),
  lbf_variable7 = numeric(),
  lbf_variable8 = numeric(),
  lbf_variable9 = numeric(),
  lbf_variable10 = numeric()
)
empty_lbf_df
}

InitEmptyInCSVariantDf <- function() {
empty_in_cs_variant_df = dplyr::tibble(
  molecular_trait_id = character(),
  variant = character(),
  chromosome = character(),
  position = integer(),
  ref = character(),
  alt = character(),
  cs_id = character(),
  cs_index = character(),
  region = character(),
  pip = numeric(),
  z = numeric(),
  cs_min_r2 = numeric(),
  cs_avg_r2 = numeric(),
  cs_size = integer(),
  posterior_mean = numeric(),
  posterior_sd = numeric(),
  cs_log10bf = numeric()
    )
empty_in_cs_variant_df
}

InitEmptyCS <- function() {
empty_cs_df = dplyr::tibble(
  molecular_trait_id = numeric(),
  cs_id = numeric(),
  cs_index = numeric(),
  region = character(),
  cs_log10bf = numeric(),
  cs_avg_r2 = numeric(),
  cs_min_r2 = numeric(),
  cs_size = numeric(),
  low_purity = numeric()
)
empty_cs_df

}
