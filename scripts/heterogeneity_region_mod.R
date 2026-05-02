# heterogeneity_region_mod.R
# Revised district-level heterogeneity analysis.
# Key change: temperature and precipitation interactions are estimated jointly.
# Equation implemented:
#   Y_ismt = beta_T T_ismt + beta_P P_ismt + rho Z_it
#            + lambda_T (T_ismt x Z_it) + lambda_P (P_ismt x Z_it)
#            + mu_im + delta_st + e_ismt
# For time-invariant Z_i, the main effect is absorbed by district-by-month FE
# and is omitted from the right-hand side.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(lfe)
  library(tibble)
})

source_candidates <- c("scripts/functions_mod.R", "functions_mod.R")
source_found <- source_candidates[file.exists(source_candidates)][1]
if (!is.na(source_found)) source(source_found)

if (!exists("weighted_center")) {
  weighted_center <- function(x, w = NULL) {
    center <- if (is.null(w)) mean(x, na.rm = TRUE) else weighted.mean(x, w, na.rm = TRUE)
    x - center
  }
}
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

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
base_path <- "inputs"
out_dir <- file.path(base_path, "outputs_heterogeneity_region_mod")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

f_base <- file.path(base_path, "SuicideData_ROK.csv")
f_grdp <- file.path(base_path, "grdp.csv")
stopifnot(file.exists(f_base))

# -----------------------------------------------------------------------------
# Load and prepare base panel
# -----------------------------------------------------------------------------
dat <- fread(f_base) %>% as.data.frame()

dat <- dat %>%
  mutate(
    fips = as.character(fips),
    year = as.integer(ifelse("year" %in% names(.), year, yr)),
    month = as.integer(month),
    state = if ("state" %in% names(.)) as.character(state) else substr(fips, 1, 2),
    fipsmo = interaction(fips, month, drop = TRUE),
    stateyear = interaction(state, year, drop = TRUE),
    ln_pop = log(pop)
  )

# -----------------------------------------------------------------------------
# Merge GRDP if available
# -----------------------------------------------------------------------------
if (file.exists(f_grdp)) {
  grdp <- fread(f_grdp) %>% as.data.frame()
  names(grdp)[1] <- "fips"
  grdp_long <- grdp %>%
    pivot_longer(cols = -fips, names_to = "year", values_to = "grdp") %>%
    mutate(
      fips = as.character(fips),
      year = as.integer(year),
      grdp = as.numeric(grdp),
      ln_grdp = log(grdp)
    ) %>%
    select(fips, year, grdp, ln_grdp)
  dat <- dat %>% left_join(grdp_long, by = c("fips", "year"))
} else {
  warning("inputs/grdp.csv was not found. Income heterogeneity will be skipped.")
  dat$ln_grdp <- NA_real_
}

# -----------------------------------------------------------------------------
# Construct district-level and centered moderators
# -----------------------------------------------------------------------------
region_means <- dat %>%
  group_by(fips) %>%
  summarise(
    mean_pop = mean(pop, na.rm = TRUE),
    avg_temp = mean(tmean, na.rm = TRUE),
    avg_prec = mean(prec, na.rm = TRUE),
    .groups = "drop"
  )

pop_threshold <- mean(dat$pop, na.rm = TRUE)
temp_threshold <- mean(dat$tmean, na.rm = TRUE)
prec_threshold <- mean(dat$prec, na.rm = TRUE)

region_means <- region_means %>%
  mutate(
    high_pop = as.integer(mean_pop >= pop_threshold),
    warm_region = as.integer(avg_temp >= temp_threshold),
    wet_region = as.integer(avg_prec >= prec_threshold)
  )

dat <- dat %>%
  left_join(region_means, by = "fips") %>%
  mutate(
    ln_grdp_c = weighted_center(ln_grdp, popw),
    ln_pop_c = weighted_center(ln_pop, popw),
    avg_temp_c = weighted_center(avg_temp, popw),
    avg_prec_c = weighted_center(avg_prec, popw)
  )

thresholds <- tibble(
  threshold = c("population", "long_term_temperature", "long_term_precipitation"),
  value = c(pop_threshold, temp_threshold, prec_threshold)
)
write_csv(thresholds, file.path(out_dir, "thresholds_used.csv"))

# -----------------------------------------------------------------------------
# Model runner
# -----------------------------------------------------------------------------
run_region_model <- function(data, moderator, label, include_main = TRUE,
                             moderator_type = c("continuous", "binary")) {
  moderator_type <- match.arg(moderator_type)
  rhs <- if (include_main) {
    paste("tmean", "prec", moderator,
          paste0("tmean:", moderator), paste0("prec:", moderator), sep = " + ")
  } else {
    paste("tmean", "prec", paste0("tmean:", moderator), paste0("prec:", moderator), sep = " + ")
  }
  fmla <- as.formula(paste0("rate_adj ~ ", rhs, " | fipsmo + stateyear | 0 | fips"))
  dsub <- data %>%
    filter(!is.na(rate_adj), !is.na(tmean), !is.na(prec), !is.na(.data[[moderator]]), !is.na(popw))
  mod <- felm(fmla, data = dsub, weights = dsub$popw)
  b <- coef(mod)
  tx <- get_term(b, c(paste0("tmean:", moderator), paste0(moderator, ":tmean")))
  px <- get_term(b, c(paste0("prec:", moderator), paste0(moderator, ":prec")))

  coef_rows <- bind_rows(
    coef_extract(mod, "tmean", "Temperature"),
    coef_extract(mod, "prec", "Precipitation"),
    coef_extract(mod, tx, "Temperature x moderator"),
    coef_extract(mod, px, "Precipitation x moderator")
  ) %>%
    mutate(moderator = label, moderator_var = moderator, N = mod$N, .before = 1)

  if (moderator_type == "binary") {
    eval_points <- tibble(group = c("Reference / below", "Above / high"), x = c(0, 1))
  } else {
    probs <- quantile(dsub[[moderator]], probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
    eval_points <- tibble(group = c("p25", "p50", "p75"), x = as.numeric(probs))
  }

  effects <- eval_points %>%
    rowwise() %>%
    mutate(
      tmp = list(lincom_effect(mod, "tmean", tx, x = x, scale = 1)),
      prc = list(lincom_effect(mod, "prec", px, x = x, scale = 100))
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
    mutate(moderator = label, moderator_var = moderator, N = mod$N, .before = 1)

  list(model = mod, coefficients = coef_rows, effects = effects)
}

specs <- tribble(
  ~moderator, ~label, ~include_main, ~type,
  "ln_grdp_c", "Log GRDP per capita", TRUE, "continuous",
  "ln_pop_c", "Log population", TRUE, "continuous",
  "high_pop", "Above-mean population", FALSE, "binary",
  "avg_temp_c", "Long-term average temperature", FALSE, "continuous",
  "avg_prec_c", "Long-term average precipitation", FALSE, "continuous",
  "warm_region", "Above-mean long-term temperature", FALSE, "binary",
  "wet_region", "Above-mean long-term precipitation", FALSE, "binary"
)

results <- vector("list", nrow(specs))
for (i in seq_len(nrow(specs))) {
  s <- specs[i, ]
  if (all(is.na(dat[[s$moderator]]))) next
  results[[i]] <- run_region_model(
    dat, moderator = s$moderator, label = s$label,
    include_main = s$include_main, moderator_type = s$type
  )
}
results <- Filter(Negate(is.null), results)

coef_tbl <- bind_rows(lapply(results, `[[`, "coefficients"))
effect_tbl <- bind_rows(lapply(results, `[[`, "effects"))

write_csv(coef_tbl, file.path(out_dir, "regional_model_coefficients.csv"))
write_csv(effect_tbl, file.path(out_dir, "regional_marginal_effects.csv"))

# Table-ready pieces: one file for temperature heterogeneity and one for precipitation heterogeneity.
table3_inputs <- coef_tbl %>%
  filter(term %in% c("Temperature", "Temperature x moderator")) %>%
  select(moderator, term, est, se, t, p, N)

table4_inputs <- coef_tbl %>%
  filter(term %in% c("Precipitation", "Precipitation x moderator")) %>%
  select(moderator, term, est, se, t, p, N)

write_csv(table3_inputs, file.path(out_dir, "table3_temperature_regional_inputs.csv"))
write_csv(table4_inputs, file.path(out_dir, "table4_precipitation_regional_inputs.csv"))

print(table3_inputs)
print(table4_inputs)
