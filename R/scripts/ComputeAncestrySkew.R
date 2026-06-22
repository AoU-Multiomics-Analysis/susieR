library(tidyverse)
library(data.table)
library(optparse)

# Compute enrichment of high-PIP variants in the GVS subpopulation with the
# highest minor allele frequency, both before and after removing admixed groups.

option_list <- list(
    optparse::make_option(c("--AnnotationData"), type = "character", default = NULL,
                        help = "Annotation data for variants containing ancestry specific AC/AN and cohort level AC/AN", metavar = "type"),
    optparse::make_option(c("--OutputPrefix"), type = "character", default = NULL,
                        help = "Output file name prefix", metavar = "type"),
    optparse::make_option(c("--PipThreshold"), type = "double", default = 0.9,
                        help = "Minimum PIP threshold used to select variants [default %default]", metavar = "number"),
    optparse::make_option(c("--AdmixedSubpops"), type = "character", default = "oth",
                        help = "Comma-separated GVS subpopulation labels to remove for the no-admixed skew calculation [default %default]", metavar = "string"),
    optparse::make_option(c("--KeepInputColumns"), action = "store_true", default = FALSE,
                        help = "Keep all input annotation columns and append ancestry skew columns [default %default]")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

if (is.null(opt$AnnotationData)) {
    stop("--AnnotationData is required", call. = FALSE)
}
if (is.null(opt$OutputPrefix)) {
    stop("--OutputPrefix is required", call. = FALSE)
}

AnnotationPath <- opt$AnnotationData
OutputName <- paste0(opt$OutputPrefix, ".AncestrySkew.tsv.gz")
PipThreshold <- opt$PipThreshold
AdmixedSubpops <- opt$AdmixedSubpops %>%
    str_split(",", simplify = FALSE) %>%
    unlist() %>%
    str_trim() %>%
    str_to_lower()
AdmixedSubpops <- AdmixedSubpops[AdmixedSubpops != ""]
KeepInputColumns <- opt$KeepInputColumns

required_columns <- c("variant", "pip", "gvs_all_ac", "gvs_all_an")

# Convert allele frequency to minor allele frequency before choosing the max
# subpopulation. This fixes older annotations where max subpop was based on AF.
maf <- function(af) {
    pmin(af, 1 - af)
}

# Report which required columns are absent from an input table.
missing_columns <- function(df, cols) {
    setdiff(cols, names(df))
}

# Detect gzip input from the magic bytes rather than relying on the file suffix.
is_gzip_file <- function(path) {
    con <- file(path, "rb")
    on.exit(close(con))
    magic <- readBin(con, what = "raw", n = 2)
    length(magic) == 2 && identical(as.integer(magic), c(0x1f, 0x8b))
}

# Read annotation tables whether Cromwell localized them with or without a gzip
# suffix.
read_annotation <- function(path) {
    if (is_gzip_file(path)) {
        return(fread(cmd = paste("gzip -dc", shQuote(path))))
    }
    fread(path)
}

# Stop with a context-specific message when expected annotation columns are
# missing.
require_columns <- function(df, cols, context) {
    missing <- missing_columns(df, cols)
    if (length(missing) > 0) {
        stop(
            paste0(context, " requires missing column(s): ", paste(missing, collapse = ", ")),
            call. = FALSE
        )
    }
}

# Translate a GVS allele-frequency column name into its subpopulation label.
get_subpop_from_af_col <- function(col) {
    col %>%
        str_remove("^gvs_") %>%
        str_remove("_af$")
}

# Select the subpopulation with the largest MAF for each variant and write both
# the chosen label and the selected MAF into prefix-specific output columns.
choose_max_subpop <- function(df, subpops, prefix) {
    maf_cols <- paste0("gvs_", subpops, "_maf")
    max_maf_col <- paste0(prefix, "_max_maf")
    max_subpop_col <- paste0(prefix, "_max_subpop")

    df_with_ids <- df %>%
        select(-any_of(c(max_subpop_col, max_maf_col))) %>%
        mutate(.row_id = row_number())

    max_subpop_by_row <- df_with_ids %>%
        select(.row_id, all_of(maf_cols)) %>%
        pivot_longer(
            cols = all_of(maf_cols),
            names_to = "maf_column",
            values_to = max_maf_col
        ) %>%
        mutate(!!max_subpop_col := get_subpop_from_af_col(str_remove(maf_column, "_maf$"))) %>%
        filter(!is.na(.data[[max_maf_col]])) %>%
        slice_max(
            order_by = .data[[max_maf_col]],
            n = 1,
            with_ties = FALSE,
            by = .row_id
        ) %>%
        select(.row_id, all_of(c(max_subpop_col, max_maf_col)))

    df_with_ids %>%
        left_join(max_subpop_by_row, by = ".row_id") %>%
        select(-.row_id)
}

# Given the selected max subpopulation, copy that subpopulation's AC/AN into
# prefix-specific max count columns.
add_max_counts <- function(df, prefix) {
    max_subpop_col <- paste0(prefix, "_max_subpop")
    max_ac_col <- paste0(prefix, "_max_ac")
    max_an_col <- paste0(prefix, "_max_an")
    df_without_existing <- df %>%
        select(-any_of(c(max_ac_col, max_an_col)))
    count_cols <- names(df_without_existing) %>% keep(~ str_detect(.x, "^gvs_[^_]+_(ac|an)$"))

    df_with_ids <- df_without_existing %>%
        mutate(.row_id = row_number()) %>%
        select(.row_id, all_of(max_subpop_col), all_of(count_cols))

    max_counts_by_row <- df_with_ids %>%
        pivot_longer(
            cols = all_of(count_cols),
            names_to = c("gvs_prefix", "subpop", ".value"),
            names_pattern = "^(gvs)_([^_]+)_(ac|an)$"
        ) %>%
        filter(subpop == .data[[max_subpop_col]]) %>%
        transmute(.row_id, !!max_ac_col := ac, !!max_an_col := an)

    df_without_existing %>%
        mutate(.row_id = row_number()) %>%
        left_join(max_counts_by_row, by = ".row_id") %>%
        select(-.row_id)
}

# Define the Fisher-test background as all cohort counts minus the max
# subpopulation counts.
add_background_counts <- function(df, prefix, all_ac_col, all_an_col) {
    max_ac_col <- paste0(prefix, "_max_ac")
    max_an_col <- paste0(prefix, "_max_an")
    bg_ac_col <- paste0(prefix, "_background_ac")
    bg_an_col <- paste0(prefix, "_background_an")

    df %>%
        mutate(
            !!bg_ac_col := .data[[all_ac_col]] - .data[[max_ac_col]],
            !!bg_an_col := .data[[all_an_col]] - .data[[max_an_col]]
        )
}

# Run one Fisher exact test per variant, returning odds ratio and p-value columns
# named for the current analysis prefix. A pseudocount of one keeps odds ratios
# finite when one valid cell is zero.
run_fisher <- function(df, prefix) {
    max_ac_col <- paste0(prefix, "_max_ac")
    max_an_col <- paste0(prefix, "_max_an")
    bg_ac_col <- paste0(prefix, "_background_ac")
    bg_an_col <- paste0(prefix, "_background_an")
    odds_col <- paste0(prefix, "_odds_ratio")
    p_col <- paste0(prefix, "_p_value")

    fisher_result <- function(max_ac, max_an, background_ac, background_an) {
        values <- c(max_ac, max_an, background_ac, background_an)

        if (any(is.na(values)) || any(values < 0)) {
            return(tibble(odds_ratio = NA_real_, p_value = NA_real_))
        }

        values <- values + 1
        test <- fisher.test(matrix(values, nrow = 2, byrow = TRUE))
        tibble(odds_ratio = unname(test$estimate), p_value = test$p.value)
    }

    df %>%
        mutate(
            .fisher = pmap(
                list(
                    .data[[max_ac_col]],
                    .data[[max_an_col]],
                    .data[[bg_ac_col]],
                    .data[[bg_an_col]]
                ),
                fisher_result
            )
        ) %>%
        unnest_wider(.fisher) %>%
        rename(!!odds_col := odds_ratio, !!p_col := p_value)
}

####### LOAD DATA ########
# Read the annotated SuSiE table and validate the minimum columns needed for
# ancestry skew calculations.
AnnotationDf <- read_annotation(AnnotationPath)

require_columns(AnnotationDf, required_columns, "Ancestry skew")

# Discover all ancestry-specific GVS AF columns present in the annotation file.
af_cols <- names(AnnotationDf) %>%
    keep(~ str_detect(.x, "^gvs_[^_]+_af$") && .x != "gvs_all_af" && .x != "gvs_max_af")

if (length(af_cols) == 0) {
    stop("No GVS subpopulation allele frequency columns found. Expected columns like gvs_afr_af.", call. = FALSE)
}

subpops <- map_chr(af_cols, get_subpop_from_af_col)
admixed_subpops <- intersect(AdmixedSubpops, subpops)
nonadmixed_subpops <- setdiff(subpops, admixed_subpops)

# Validate that the requested admixed labels exist and that at least one
# comparison population remains after removing them.
if (length(admixed_subpops) == 0) {
    stop(
        paste0(
            "None of the requested admixed subpops were found in GVS AF columns: ",
            paste(AdmixedSubpops, collapse = ", ")
        ),
        call. = FALSE
    )
}
if (length(nonadmixed_subpops) == 0) {
    stop("No non-admixed subpopulations remain after applying --AdmixedSubpops.", call. = FALSE)
}

count_cols <- c(paste0("gvs_", subpops, "_ac"), paste0("gvs_", subpops, "_an"))
require_columns(
    AnnotationDf,
    count_cols,
    paste0(
        "MAF-based max subpopulation AC/AN lookup. The attached annotation contains AF columns, ",
        "but exact skew needs matching per-subpopulation AC/AN columns"
    )
)

numeric_cols <- names(AnnotationDf) %>%
    keep(~ str_detect(.x, "^gvs_.*_(ac|an|af)$") || .x %in% c("pip"))
numeric_cols <- unique(numeric_cols)

# Keep output columns explicit so empty and non-empty outputs have the same
# schema. These are the only columns added by ancestry-skew mode when
# --KeepInputColumns is used.
OutputColumns <- c(
    "variant", "pip",
    "gvs_max_subpop", "gvs_max_af",
    "gvs_odds_ratio", "gvs_p_value",
    "gvs_no_admixed_odds_ratio", "gvs_no_admixed_p_value"
)

# Keep only variants at or above the requested PIP threshold, then create MAF
# columns for every GVS ancestry-specific AF column.
SkewInput <- AnnotationDf %>%
    mutate(across(all_of(numeric_cols), ~ as.numeric(.))) %>%
    filter(pip >= PipThreshold) %>%
    mutate(across(all_of(af_cols), maf, .names = "{.col}_maf"))

maf_col_names <- paste0(af_cols, "_maf")
names(SkewInput)[match(maf_col_names, names(SkewInput))] <- paste0("gvs_", subpops, "_maf")

if (nrow(SkewInput) == 0) {
    warning("No variants passed the PIP threshold; writing an empty output.")
    empty_columns <- OutputColumns
    if (KeepInputColumns) {
        empty_columns <- c(names(AnnotationDf), setdiff(OutputColumns, names(AnnotationDf)))
    }
    empty_output <- as_tibble(setNames(replicate(length(empty_columns), logical(0), simplify = FALSE), empty_columns))
    empty_output %>% write_tsv(OutputName)
    quit(save = "no", status = 0)
}

# Round 1: ancestry skew with admixed samples included.
#
# This recomputes the max GVS subpopulation from MAF, looks up that
# subpopulation's AC/AN, builds the all-other-background counts, and runs the
# first Fisher test.
SkewInput <- SkewInput %>%
    choose_max_subpop(subpops, "gvs") %>%
    mutate(gvs_max_af = gvs_max_maf) %>%
    add_max_counts("gvs") %>%
    add_background_counts("gvs", "gvs_all_ac", "gvs_all_an") %>%
    run_fisher("gvs")

# Round 2: ancestry skew after admixed samples are removed.
#
# The admixed AC/AN are subtracted from cohort-level AC/AN first. Then the max
# subpopulation is recalculated from MAF using only the non-admixed
# subpopulations before running the second Fisher test.
SkewInput <- SkewInput %>%
    mutate(
        gvs_no_admixed_all_ac = gvs_all_ac - rowSums(across(all_of(paste0("gvs_", admixed_subpops, "_ac"))), na.rm = FALSE),
        gvs_no_admixed_all_an = gvs_all_an - rowSums(across(all_of(paste0("gvs_", admixed_subpops, "_an"))), na.rm = FALSE)
    ) %>%
    choose_max_subpop(nonadmixed_subpops, "gvs_no_admixed") %>%
    add_max_counts("gvs_no_admixed") %>%
    add_background_counts("gvs_no_admixed", "gvs_no_admixed_all_ac", "gvs_no_admixed_all_an") %>%
    run_fisher("gvs_no_admixed")

OutputDf <- SkewInput %>%
    select(all_of(OutputColumns))

# Optionally preserve the full annotation input beside the newly computed skew
# columns for downstream ad hoc review.
if (KeepInputColumns) {
    OutputDf <- SkewInput %>%
        select(all_of(names(AnnotationDf)), all_of(setdiff(OutputColumns, names(AnnotationDf))))
}

OutputDf %>%
    write_tsv(OutputName)
