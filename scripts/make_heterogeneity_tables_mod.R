# make_heterogeneity_tables_mod.R
# -----------------------------------------------------------------------------
# Purpose
#   Build publication-ready CSV/Excel tables from the revised heterogeneity
#   analysis outputs.
#
# Required upstream scripts, run before this file:
#   1. heterogeneity_region_mod.R
#   2. sex_mod.R
#   3. age_mod.R
#
# Main outputs:
#   inputs/outputs_heterogeneity_tables_mod/Table_3_Heterogeneity_Coefficients.csv
#   inputs/outputs_heterogeneity_tables_mod/Table_4_Marginal_Effects.csv
#   inputs/outputs_heterogeneity_tables_mod/heterogeneity_tables.xlsx, if writexl
#   or openxlsx is installed.
#
# Design principle
#   Table 3 reports whether heterogeneity exists: interaction coefficients.
#   Table 4 reports what the estimates imply: marginal effects at intuitive
#   moderator values or demographic groups.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(tibble)
})

base_path <- "inputs"
out_dir <- file.path(base_path, "outputs_heterogeneity_tables_mod")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
paths <- list(
  region_coef = file.path(base_path, "outputs_heterogeneity_region_mod", "regional_model_coefficients.csv"),
  region_eff  = file.path(base_path, "outputs_heterogeneity_region_mod", "regional_marginal_effects.csv"),
  sex_coef    = file.path(base_path, "outputs_sex_mod", "sex_model_coefficients.csv"),
  sex_eff     = file.path(base_path, "outputs_sex_mod", "sex_marginal_effects.csv"),
  age_coef    = file.path(base_path, "outputs_age_mod", "age_model_coefficients.csv"),
  age_eff_rep = file.path(base_path, "outputs_age_mod", "representative_age_group_effects.csv"),
  age_eff_all = file.path(base_path, "outputs_age_mod", "age_marginal_effects_by_agebin.csv")
)

required <- paths[c("region_coef", "region_eff", "sex_coef", "sex_eff", "age_coef", "age_eff_rep")]
missing_required <- names(required)[!file.exists(unlist(required))]
if (length(missing_required) > 0) {
  stop(
    "Missing upstream output files: ", paste(missing_required, collapse = ", "),
    "\nRun heterogeneity_region_mod.R, sex_mod.R, and age_mod.R before this script."
  )
}

read_csv_quiet <- function(path) {
  readr::read_csv(path, show_col_types = FALSE)
}

region_coef <- read_csv_quiet(paths$region_coef)
region_eff  <- read_csv_quiet(paths$region_eff)
sex_coef    <- read_csv_quiet(paths$sex_coef)
sex_eff     <- read_csv_quiet(paths$sex_eff)
age_coef    <- read_csv_quiet(paths$age_coef)
age_eff_rep <- read_csv_quiet(paths$age_eff_rep)
age_eff_all <- if (file.exists(paths$age_eff_all)) read_csv_quiet(paths$age_eff_all) else NULL

# -----------------------------------------------------------------------------
# Formatting helpers
# -----------------------------------------------------------------------------
star <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "",
    p < 0.01 ~ "***",
    p < 0.05 ~ "**",
    p < 0.10 ~ "*",
    TRUE ~ ""
  )
}

fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), "", formatC(x, format = "f", digits = digits))
}

fmt_est_se <- function(est, se, p, digits = 3) {
  ifelse(
    is.na(est),
    "",
    paste0(fmt_num(est, digits), star(p), " (", fmt_num(se, digits), ")")
  )
}

fmt_ci <- function(est, lower, upper, digits = 3) {
  ifelse(
    is.na(est),
    "",
    paste0(fmt_num(est, digits), " [", fmt_num(lower, digits), ", ", fmt_num(upper, digits), "]")
  )
}

fmt_pct <- function(x, digits = 2) {
  ifelse(is.na(x), "", paste0(fmt_num(x, digits), "%"))
}

row_term <- function(df, term_name, scale = 1) {
  z <- df %>% filter(.data$term == term_name) %>% slice(1)
  if (nrow(z) == 0) {
    return(tibble(est = NA_real_, se = NA_real_, p = NA_real_, N = NA_real_))
  }
  tibble(
    est = z$est * scale,
    se = z$se * abs(scale),
    p = z$p,
    N = z$N
  )
}

make_coef_row <- function(df, moderator_label, panel, equation, reference_note, fe_note, weight_note) {
  temp_base <- row_term(df, "Temperature", scale = 1)
  temp_int <- df %>% filter(str_detect(.data$term, regex("Temperature x", ignore_case = TRUE))) %>% slice(1)
  precip_base <- row_term(df, "Precipitation", scale = 100)
  precip_int <- df %>% filter(str_detect(.data$term, regex("Precipitation x", ignore_case = TRUE))) %>% slice(1)

  if (nrow(temp_int) == 0) temp_int <- tibble(est = NA_real_, se = NA_real_, p = NA_real_, N = NA_real_)
  if (nrow(precip_int) == 0) precip_int <- tibble(est = NA_real_, se = NA_real_, p = NA_real_, N = NA_real_)

  # Scale precipitation interaction to a 100 mm change for readability.
  precip_int <- precip_int %>% mutate(est = .data$est * 100, se = .data$se * 100)

  N_use <- suppressWarnings(max(c(temp_base$N, precip_base$N, temp_int$N, precip_int$N), na.rm = TRUE))
  if (!is.finite(N_use)) N_use <- NA_real_

  tibble(
    panel = panel,
    equation = equation,
    moderator = moderator_label,
    reference_or_center = reference_note,
    temperature_effect_at_reference = fmt_est_se(temp_base$est, temp_base$se, temp_base$p, digits = 4),
    temperature_x_moderator = fmt_est_se(temp_int$est, temp_int$se, temp_int$p, digits = 4),
    precipitation_effect_at_reference_per_100mm = fmt_est_se(precip_base$est, precip_base$se, precip_base$p, digits = 4),
    precipitation_x_moderator_per_100mm = fmt_est_se(precip_int$est, precip_int$se, precip_int$p, digits = 4),
    observations = N_use,
    fixed_effects = fe_note,
    weights = weight_note
  )
}

# -----------------------------------------------------------------------------
# Table 3. Interaction coefficients
# -----------------------------------------------------------------------------
region_table3 <- region_coef %>%
  group_by(.data$moderator) %>%
  group_split() %>%
  lapply(function(df) {
    lab <- unique(df$moderator)[1]
    ref_note <- if (any(str_detect(df$moderator_var, "_c$"))) {
      "Continuous moderator centered at its weighted mean"
    } else if (any(df$moderator_var %in% c("high_pop", "warm_region", "wet_region"))) {
      "Reference group is below-threshold district"
    } else {
      "See upstream script"
    }
    make_coef_row(
      df = df,
      moderator_label = lab,
      panel = "A. District-level heterogeneity",
      equation = "Equation (2)",
      reference_note = ref_note,
      fe_note = "District x month; upper-region x year",
      weight_note = "District population weight (popw)"
    )
  }) %>%
  bind_rows()

sex_table3 <- make_coef_row(
  df = sex_coef,
  moderator_label = "Male (reference: female)",
  panel = "B. Demographic heterogeneity",
  equation = "Equation (3)",
  reference_note = "Female is the reference group",
  fe_note = "District x month x sex; upper-region x year x sex",
  weight_note = "Sex-specific population"
)

age_table3 <- make_coef_row(
  df = age_coef,
  moderator_label = "Age midpoint, centered and divided by 5",
  panel = "B. Demographic heterogeneity",
  equation = "Equation (3)",
  reference_note = "One moderator unit equals a 5-year increase from the weighted mean age midpoint",
  fe_note = "District x month x age group; upper-region x year x age group",
  weight_note = "Age-group-specific population"
)

table3 <- bind_rows(region_table3, sex_table3, age_table3)

# -----------------------------------------------------------------------------
# Table 4. Implied marginal effects
# -----------------------------------------------------------------------------
region_table4 <- region_eff %>%
  transmute(
    panel = "A. District-level heterogeneity",
    moderator = .data$moderator,
    evaluation_point = dplyr::recode(.data$group, p25 = "25th percentile", p50 = "Median", p75 = "75th percentile", .default = .data$group),
    moderator_value_used = .data$x,
    baseline_rate = NA_real_,
    temperature_effect_per_1C = fmt_ci(.data$temp_effect_per_1C, .data$temp_lower, .data$temp_upper, digits = 4),
    temperature_percent_change = "",
    precipitation_effect_per_100mm = fmt_ci(.data$precip_effect_per_100mm, .data$precip_lower, .data$precip_upper, digits = 4),
    precipitation_percent_change = "",
    observations = .data$N,
    notes = "Continuous moderator values are centered values; binary moderator values are 0/1."
  )

sex_table4 <- sex_eff %>%
  transmute(
    panel = "B. Demographic heterogeneity",
    moderator = "Sex",
    evaluation_point = .data$sex_label,
    moderator_value_used = .data$male,
    baseline_rate = .data$baseline_rate,
    temperature_effect_per_1C = fmt_ci(.data$temp_effect_per_1C, .data$temp_lower, .data$temp_upper, digits = 4),
    temperature_percent_change = fmt_pct(.data$temp_pct_per_1C, digits = 2),
    precipitation_effect_per_100mm = fmt_ci(.data$precip_effect_per_100mm, .data$precip_lower, .data$precip_upper, digits = 4),
    precipitation_percent_change = fmt_pct(.data$precip_pct_per_100mm, digits = 2),
    observations = NA_real_,
    notes = "Effects are evaluated for female and male suicide rates separately."
  )

age_table4 <- age_eff_rep %>%
  transmute(
    panel = "B. Demographic heterogeneity",
    moderator = "Age",
    evaluation_point = .data$age_group,
    moderator_value_used = .data$age5_c,
    baseline_rate = .data$baseline_rate,
    temperature_effect_per_1C = fmt_ci(.data$temp_effect_per_1C, .data$temp_lower, .data$temp_upper, digits = 4),
    temperature_percent_change = fmt_pct(.data$temp_pct_per_1C, digits = 2),
    precipitation_effect_per_100mm = fmt_ci(.data$precip_effect_per_100mm, .data$precip_lower, .data$precip_upper, digits = 4),
    precipitation_percent_change = fmt_pct(.data$precip_pct_per_100mm, digits = 2),
    observations = NA_real_,
    notes = "Representative age groups are evaluated at population-weighted mean age5_c within each group."
  )

table4 <- bind_rows(region_table4, sex_table4, age_table4)

# Optional appendix: all age-bin effects.
if (!is.null(age_eff_all)) {
  appendix_age <- age_eff_all %>%
    transmute(
      age_code = .data$age_code,
      age_label = .data$age_label,
      age_mid = .data$age_mid,
      baseline_rate = .data$baseline_rate,
      temperature_effect_per_1C = fmt_ci(.data$temp_effect_per_1C, .data$temp_lower, .data$temp_upper, digits = 4),
      temperature_percent_change = fmt_pct(.data$temp_pct_per_1C, digits = 2),
      precipitation_effect_per_100mm = fmt_ci(.data$precip_effect_per_100mm, .data$precip_lower, .data$precip_upper, digits = 4),
      precipitation_percent_change = fmt_pct(.data$precip_pct_per_100mm, digits = 2)
    )
} else {
  appendix_age <- tibble()
}

# Also save unformatted numeric versions for replication and downstream editing.
raw_coefficients <- bind_rows(
  region_coef %>% mutate(source = "district_level"),
  sex_coef %>% mutate(source = "sex"),
  age_coef %>% mutate(source = "age")
)

raw_marginal_effects <- bind_rows(
  region_eff %>% mutate(source = "district_level"),
  sex_eff %>% mutate(source = "sex"),
  age_eff_rep %>% mutate(source = "age_representative_groups")
)

# README sheet/table.
readme_tbl <- tibble(
  item = c(
    "Table 3 purpose",
    "Table 4 purpose",
    "Temperature scaling",
    "Precipitation scaling",
    "Significance stars",
    "District-level FE",
    "Demographic FE",
    "CSV outputs",
    "Excel output"
  ),
  description = c(
    "Interaction-coefficient table: shows whether weather effects differ by moderator.",
    "Marginal-effect table: shows implied effects at selected moderator values or demographic groups.",
    "Temperature effects are reported per 1 degree C increase in monthly mean temperature.",
    "Precipitation effects are reported per 100 mm increase in monthly precipitation for readability, although models use precipitation in mm.",
    "*** p < 0.01, ** p < 0.05, * p < 0.10.",
    "District-level models use district-by-month and upper-region-by-year fixed effects.",
    "Sex and age models use group-specific district-by-month and upper-region-by-year fixed effects.",
    "This script always writes CSV files.",
    "This script writes an XLSX workbook if package writexl or openxlsx is installed."
  )
)

# -----------------------------------------------------------------------------
# Write outputs
# -----------------------------------------------------------------------------
write_csv(table3, file.path(out_dir, "Table_3_Heterogeneity_Coefficients.csv"))
write_csv(table4, file.path(out_dir, "Table_4_Marginal_Effects.csv"))
write_csv(raw_coefficients, file.path(out_dir, "Raw_Heterogeneity_Coefficients.csv"))
write_csv(raw_marginal_effects, file.path(out_dir, "Raw_Marginal_Effects.csv"))
write_csv(readme_tbl, file.path(out_dir, "README_for_tables.csv"))
if (nrow(appendix_age) > 0) {
  write_csv(appendix_age, file.path(out_dir, "Appendix_AgeBin_Marginal_Effects.csv"))
}

xlsx_path <- file.path(out_dir, "heterogeneity_tables.xlsx")
excel_sheets <- list(
  README = readme_tbl,
  Table3_Coefficients = table3,
  Table4_MarginalEffects = table4,
  Raw_Coefficients = raw_coefficients,
  Raw_MarginalEffects = raw_marginal_effects
)
if (nrow(appendix_age) > 0) {
  excel_sheets$Appendix_AgeBins <- appendix_age
}

if (requireNamespace("writexl", quietly = TRUE)) {
  writexl::write_xlsx(excel_sheets, path = xlsx_path)
  message("Wrote Excel workbook: ", xlsx_path)
} else if (requireNamespace("openxlsx", quietly = TRUE)) {
  openxlsx::write.xlsx(excel_sheets, file = xlsx_path, overwrite = TRUE)
  message("Wrote Excel workbook: ", xlsx_path)
} else {
  message("Package writexl/openxlsx not installed. CSV files were written, but XLSX was skipped.")
  message("To enable XLSX output, run install.packages('writexl') and rerun this script.")
}

message("Done. Table outputs are in: ", out_dir)
message("Main files:")
message("  - ", file.path(out_dir, "Table_3_Heterogeneity_Coefficients.csv"))
message("  - ", file.path(out_dir, "Table_4_Marginal_Effects.csv"))
