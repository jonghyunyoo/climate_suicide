suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(lfe)
})

# -----------------------------
# 1) suicide data
# -----------------------------
dat <- fread("inputs/SuicideData_ROK.csv") %>% as.data.frame()

dat <- dat %>%
  mutate(
    fips = as.numeric(fips),
    year = as.integer(year),
    month = as.integer(month)
  ) %>%
  arrange(fips, year, month) %>%
  mutate(
    fipsmo = fips * 100 + month,
    stateyear = floor(fips / 1000) * 10000 + year,
    ln_pop = log(pop)
  )

# -----------------------------
# 2) continuous interaction
# -----------------------------
mod_pop <- felm(
  rate_adj ~ tmean + ln_pop + tmean:ln_pop + prec |
    fipsmo + stateyear |
    0 |
    fips,
  data = dat,
  weights = dat$popw
)

summary(mod_pop)

# -----------------------------
# 3) mean population by region
# -----------------------------
pop_mean_by_fips <- dat %>%
  group_by(fips) %>%
  summarise(mean_pop = mean(pop, na.rm = TRUE), .groups = "drop")

pop_avg <- 201861

pop_mean_by_fips <- pop_mean_by_fips %>%
  mutate(
    high_pop = ifelse(mean_pop >= pop_avg, 1, 0)
  )

# -----------------------------
# 4) merge group dummy
# -----------------------------
dat_avg <- dat %>%
  left_join(pop_mean_by_fips %>% select(fips, high_pop), by = "fips")

# -----------------------------
# 5) grouped interaction
# -----------------------------
mod_pop_avg <- felm(
  rate_adj ~ tmean + high_pop + tmean:high_pop + prec |
    fipsmo + stateyear |
    0 |
    fips,
  data = dat_avg,
  weights = dat_avg$popw
)

summary(mod_pop_avg)

# -----------------------------
# 6) group-specific effects
# -----------------------------
b <- coef(mod_pop_avg)
V <- vcov(mod_pop_avg)

int_name <- if ("tmean:high_pop" %in% names(b)) {
  "tmean:high_pop"
} else if ("high_pop:tmean" %in% names(b)) {
  "high_pop:tmean"
} else {
  stop("상호작용항 이름을 찾을 수 없습니다.")
}

b1 <- b["tmean"]
b3 <- b[int_name]

# Low population
est_low <- b1
se_low  <- sqrt(V["tmean", "tmean"])

# High population
est_high <- b1 + b3
var_high <- V["tmean", "tmean"] + V[int_name, int_name] +
  2 * V["tmean", int_name]
se_high <- sqrt(var_high)

res_avg <- data.frame(
  group = c("Below mean population", "Above mean population"),
  estimate = c(est_low, est_high),
  se = c(se_low, se_high)
) %>%
  mutate(
    ci_low = estimate - 1.96 * se,
    ci_high = estimate + 1.96 * se
  )

print(res_avg)

baseline_df <- dat_avg %>%
  group_by(high_pop) %>%
  summarise(
    baseline = mean(rate_adj, na.rm = TRUE),
    .groups = "drop"
  )

print(baseline_df)
baseline_low  <- baseline_df$baseline[baseline_df$high_pop == 0]
baseline_high <- baseline_df$baseline[baseline_df$high_pop == 1]

pct_low  <- (est_low  / baseline_low)  * 100
pct_high <- (est_high / baseline_high) * 100
res_avg <- res_avg %>%
  mutate(
    pct_change = c(pct_low, pct_high)
  )

print(res_avg)
