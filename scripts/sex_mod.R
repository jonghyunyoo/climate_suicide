# sex_mod.R
# Revised sex-specific heterogeneity analysis.
# Key changes relative to the original script:
#   1) Temperature and precipitation interactions are estimated jointly.
#   2) The stacked sex panel uses sex-specific fixed effects.
#   3) The model is weighted by sex-specific population.
# Equation implemented:
#   Y_igsmt = beta_T T_ismt + beta_P P_ismt
#             + lambda_T (T_ismt x Male_g) + lambda_P (P_ismt x Male_g)
#             + mu_img + delta_stg + e_igsmt

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(lfe)
  library(tibble)
})

source_candidates <- c("scripts/functions_mod.R", "functions_mod.R")
source_found <- source_candidates[file.exists(source_candidates)][1]
if (!is.na(source_found)) source(source_found)

if (!exists("get_term")) {
  get_term <- function(coefs, candidates) {
    hit <- candidates[candidates %in% names(coefs)]
    if (length(hit) == 0) return(NA_character_)
    hit[1]
  }
}
if (!exists("lincom_effect")) {
  lincom_effect <- function(model, base_term, interaction_term = NULL, x = 0, scale = 1) {
    b <- coef(model); V <- vcov(model)
    est <- unname(b[base_term]); var <- V[base_term, base_term]
    if (!is.null(interaction_term) && !is.na(interaction_term) && interaction_term %in% names(b)) {
      est <- est + x * unname(b[interaction_term])
      var <- var + x^2 * V[interaction_term, interaction_term] + 2 * x * V[base_term, interaction_term]
    }
    est <- est * scale; se <- sqrt(var) * abs(scale)
    data.frame(est = est, se = se, lower = est - 1.96 * se, upper = est + 1.96 * se)
  }
}
if (!exists("coef_extract")) {
  coef_extract <- function(model, term, label = term, scale = 1) {
    b <- coef(model); V <- vcov(model)
    if (!term %in% names(b)) return(data.frame(term = label, est = NA_real_, se = NA_real_, t = NA_real_, p = NA_real_))
    se <- sqrt(V[term, term]); est <- unname(b[term]); tval <- est / se
    data.frame(term = label, est = est * scale, se = se * abs(scale), t = tval,
               p = 2 * pnorm(abs(tval), lower.tail = FALSE))
  }
}

base_path <- "inputs"
out_dir <- file.path(base_path, "outputs_sex_mod")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

f_base <- file.path(base_path, "SuicideData_ROK.csv")
f_sex <- file.path(base_path, "suicide_sex.csv")
f_pop <- file.path(base_path, "age_population.csv")
stopifnot(file.exists(f_base), file.exists(f_sex), file.exists(f_pop))

base <- read_csv(f_base, show_col_types = FALSE) %>%
  mutate(
    fips = as.character(fips),
    year = as.integer(ifelse("year" %in% names(.), year, yr)),
    month = as.integer(month),
    state = if ("state" %in% names(.)) as.character(state) else substr(fips, 1, 2),
    fipsmo = interaction(fips, month, drop = TRUE),
    stateyear = interaction(state, year, drop = TRUE)
  )

sexD <- read_csv(f_sex, show_col_types = FALSE) %>%
  mutate(fips = as.character(fips), year = as.integer(year), month = as.integer(month))

popD <- read_csv(f_pop, show_col_types = FALSE) %>%
  mutate(fips = as.character(fips), sex = as.integer(sex), age = as.integer(age))

base2 <- base %>%
  select(fips, year, month, state, fipsmo, stateyear, tmean, prec)

sex_long <- sexD %>%
  select(fips, year, month, sex_1, sex_2) %>%
  pivot_longer(cols = c(sex_1, sex_2), names_to = "sex_var", values_to = "deaths") %>%
  mutate(
    sex = ifelse(sex_var == "sex_1", 1L, 2L),
    sex_label = ifelse(sex == 1L, "Male", "Female"),
    male = ifelse(sex == 1L, 1L, 0L)
  ) %>%
  select(fips, year, month, sex, sex_label, male, deaths)

pop_long <- popD %>%
  filter(sex %in% c(1L, 2L)) %>%
  pivot_longer(cols = starts_with("pop_"), names_to = "year_var", values_to = "population") %>%
  mutate(year = as.integer(str_remove(year_var, "pop_"))) %>%
  group_by(fips, sex, year) %>%
  summarise(population = sum(population, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    sex_label = ifelse(sex == 1L, "Male", "Female"),
    male = ifelse(sex == 1L, 1L, 0L)
  ) %>%
  tidyr::crossing(month = 1:12) %>%
  mutate(month = as.integer(month))

dat <- sex_long %>%
  left_join(pop_long, by = c("fips", "year", "month", "sex", "sex_label", "male")) %>%
  left_join(base2, by = c("fips", "year", "month")) %>%
  mutate(
    deaths = as.numeric(deaths),
    population = as.numeric(population),
    rate_adj = ifelse(!is.na(population) & population > 0, deaths / population * 100000, NA_real_),
    fipsmo_sex = interaction(fipsmo, sex, drop = TRUE),
    stateyear_sex = interaction(stateyear, sex, drop = TRUE)
  ) %>%
  filter(!is.na(rate_adj), !is.na(tmean), !is.na(prec), !is.na(population), population > 0)

baseline_tbl <- dat %>%
  group_by(sex_label, male) %>%
  summarise(
    total_deaths = sum(deaths, na.rm = TRUE),
    total_pop = sum(population, na.rm = TRUE),
    baseline_rate = total_deaths / total_pop * 100000,
    mean_rate_cell = mean(rate_adj, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(baseline_tbl, file.path(out_dir, "baseline_by_sex.csv"))

mod_sex <- felm(
  rate_adj ~ tmean + prec + tmean:male + prec:male |
    fipsmo_sex + stateyear_sex | 0 | fips,
  data = dat,
  weights = dat$population
)

b <- coef(mod_sex)
tx <- get_term(b, c("tmean:male", "male:tmean"))
px <- get_term(b, c("prec:male", "male:prec"))

coef_tbl <- bind_rows(
  coef_extract(mod_sex, "tmean", "Temperature"),
  coef_extract(mod_sex, "prec", "Precipitation"),
  coef_extract(mod_sex, tx, "Temperature x Male"),
  coef_extract(mod_sex, px, "Precipitation x Male")
) %>%
  mutate(N = mod_sex$N, .before = 1)
write_csv(coef_tbl, file.path(out_dir, "sex_model_coefficients.csv"))

effect_tbl <- tibble(sex_label = c("Female", "Male"), male = c(0, 1)) %>%
  rowwise() %>%
  mutate(
    tmp = list(lincom_effect(mod_sex, "tmean", tx, x = male, scale = 1)),
    prc = list(lincom_effect(mod_sex, "prec", px, x = male, scale = 100))
  ) %>%
  ungroup() %>%
  mutate(
    temp_effect_per_1C = sapply(tmp, function(z) z$est),
    temp_se = sapply(tmp, function(z) z$se),
    temp_lower = sapply(tmp, function(z) z$lower),
    temp_upper = sapply(tmp, function(z) z$upper),
    precip_effect_per_100mm = sapply(prc, function(z) z$est),
    precip_se = sapply(prc, function(z) z$se),
    precip_lower = sapply(prc, function(z) z$lower),
    precip_upper = sapply(prc, function(z) z$upper)
  ) %>%
  select(-tmp, -prc) %>%
  left_join(baseline_tbl %>% select(sex_label, baseline_rate), by = "sex_label") %>%
  mutate(
    temp_pct_per_1C = temp_effect_per_1C / baseline_rate * 100,
    precip_pct_per_100mm = precip_effect_per_100mm / baseline_rate * 100
  )
write_csv(effect_tbl, file.path(out_dir, "sex_marginal_effects.csv"))

print(coef_tbl)
print(effect_tbl)
