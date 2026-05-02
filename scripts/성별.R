suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(lfe)
  library(ggplot2)
})

# =========================================================
# 0) 경로 설정
# =========================================================
base_path <- "inputs"
out_dir   <- file.path(base_path, "outputs_sex")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

f_base <- file.path(base_path, "SuicideData_ROK.csv")
f_sex  <- file.path(base_path, "suicide_sex.csv")
f_pop  <- file.path(base_path, "age_population.csv")

stopifnot(file.exists(f_base), file.exists(f_sex), file.exists(f_pop))

# =========================================================
# 1) 데이터 불러오기
# =========================================================
base <- readr::read_csv(f_base, show_col_types = FALSE)
sexD <- readr::read_csv(f_sex,  show_col_types = FALSE)
popD <- readr::read_csv(f_pop,  show_col_types = FALSE)

# =========================================================
# 2) 공통 키 정리
# =========================================================
base <- base %>%
  dplyr::mutate(
    fips  = as.character(fips),
    year  = as.integer(year),
    month = as.integer(month)
  )

sexD <- sexD %>%
  dplyr::mutate(
    fips  = as.character(fips),
    year  = as.integer(year),
    month = as.integer(month)
  )

popD <- popD %>%
  dplyr::mutate(
    fips = as.character(fips),
    sex  = as.integer(sex),
    age  = as.integer(age)
  )

# =========================================================
# 3) 기후/통제변수 준비
# =========================================================
if (!"state" %in% names(base)) {
  base <- base %>%
    dplyr::mutate(state = substr(fips, 1, 2))
}

base2 <- base %>%
  dplyr::mutate(
    state     = as.character(state),
    fipsmo    = interaction(fips, month, drop = TRUE),
    stateyear = interaction(state, year, drop = TRUE)
  ) %>%
  dplyr::select(fips, year, month, state, fipsmo, stateyear, tmean, prec)

# =========================================================
# 4) suicide_sex long 변환
#    sex_1 = 남자, sex_2 = 여자
# =========================================================
sex_long <- sexD %>%
  dplyr::select(fips, year, month, sex_1, sex_2) %>%
  tidyr::pivot_longer(
    cols = c(sex_1, sex_2),
    names_to = "sex_var",
    values_to = "deaths"
  ) %>%
  dplyr::mutate(
    sex       = ifelse(sex_var == "sex_1", 1L, 2L),
    sex_label = ifelse(sex == 1L, "Male", "Female"),
    male      = ifelse(sex == 1L, 1L, 0L)
  ) %>%
  dplyr::select(fips, year, month, sex, sex_label, male, deaths)

# =========================================================
# 5) age_population 처리
#    구조:
#    fips, sex, age, pop_1997 ... pop_2024
#    -> long(year, population)
#    -> sex별 총인구(연령합)
#    -> 월단위 확장
# =========================================================
pop_long <- popD %>%
  tidyr::pivot_longer(
    cols = dplyr::starts_with("pop_"),
    names_to = "year_var",
    values_to = "population"
  ) %>%
  dplyr::mutate(
    year = as.integer(stringr::str_remove(year_var, "pop_"))
  ) %>%
  dplyr::group_by(fips, sex, year) %>%
  dplyr::summarise(
    population = sum(population, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    sex_label = ifelse(sex == 1L, "Male", "Female"),
    male      = ifelse(sex == 1L, 1L, 0L)
  ) %>%
  tidyr::crossing(month = 1:12) %>%
  dplyr::mutate(month = as.integer(month)) %>%
  dplyr::select(fips, year, month, sex, sex_label, male, population)

# =========================================================
# 6) 병합 및 자살률 계산
# =========================================================
dat <- sex_long %>%
  dplyr::left_join(
    pop_long,
    by = c("fips", "year", "month", "sex", "sex_label", "male")
  ) %>%
  dplyr::left_join(
    base2,
    by = c("fips", "year", "month")
  ) %>%
  dplyr::mutate(
    deaths     = as.numeric(deaths),
    population = as.numeric(population),
    rate_adj   = ifelse(
      !is.na(population) & population > 0,
      deaths / population * 100000,
      NA_real_
    )
  ) %>%
  dplyr::filter(!is.na(rate_adj), !is.na(tmean), !is.na(prec))

cat("Final analytic sample n =", nrow(dat), "\n")
cat("Year range:", min(dat$year, na.rm = TRUE), "~", max(dat$year, na.rm = TRUE), "\n")

# =========================================================
# 7) baseline 자살률
#    정확한 baseline = 전체 deaths / 전체 population * 100,000
# =========================================================
baseline_tbl <- dat %>%
  dplyr::group_by(sex_label) %>%
  dplyr::summarise(
    total_deaths   = sum(deaths, na.rm = TRUE),
    total_pop      = sum(population, na.rm = TRUE),
    baseline_rate  = total_deaths / total_pop * 100000,
    mean_rate_cell = mean(rate_adj, na.rm = TRUE),
    .groups = "drop"
  )

print(baseline_tbl)
readr::write_csv(baseline_tbl, file.path(out_dir, "baseline_by_sex.csv"))

# =========================================================
# 8) 성별 상호작용 회귀
#    female 기준집단
# =========================================================
mod_sex <- lfe::felm(
  rate_adj ~ tmean + male + tmean:male + prec |
    fipsmo + stateyear | 0 | fips,
  data = dat
)

mod_sum <- summary(mod_sex)
print(mod_sum)

# 계수표 저장
coef_df <- as.data.frame(mod_sum$coefficients)
coef_df$term <- rownames(coef_df)
rownames(coef_df) <- NULL
readr::write_csv(coef_df, file.path(out_dir, "sex_interaction_coefficients.csv"))

# =========================================================
# 9) 성별별 1℃ 효과 계산
# =========================================================
b <- coef(mod_sex)
V <- vcov(mod_sex)

int_name <- if ("tmean:male" %in% names(b)) "tmean:male" else "male:tmean"

female_est <- unname(b["tmean"])
female_se  <- sqrt(V["tmean", "tmean"])

male_est <- unname(b["tmean"] + b[int_name])
male_var <- V["tmean", "tmean"] + V[int_name, int_name] + 2 * V["tmean", int_name]
male_se  <- sqrt(male_var)

diff_est <- unname(b[int_name])
diff_se  <- sqrt(V[int_name, int_name])

effect_tbl <- tibble::tibble(
  sex_label = c("Female", "Male", "Male - Female"),
  est = c(female_est, male_est, diff_est),
  se  = c(female_se,  male_se,  diff_se)
) %>%
  dplyr::mutate(
    lower = est - 1.96 * se,
    upper = est + 1.96 * se
  )

print(effect_tbl)
readr::write_csv(effect_tbl, file.path(out_dir, "temperature_effect_by_sex.csv"))

# =========================================================
# 10) baseline 대비 % 변화
# =========================================================
plot_tbl <- effect_tbl %>%
  dplyr::filter(sex_label %in% c("Female", "Male")) %>%
  dplyr::left_join(
    baseline_tbl %>% dplyr::select(sex_label, baseline_rate),
    by = "sex_label"
  ) %>%
  dplyr::mutate(
    pct_change = est / baseline_rate * 100,
    pct_lower  = lower / baseline_rate * 100,
    pct_upper  = upper / baseline_rate * 100
  )

print(plot_tbl)

# =========================================================
# 11) 요약표
# =========================================================
summary_tbl <- plot_tbl %>%
  dplyr::mutate(
    effect_per_1C = sprintf("%.4f (%.4f, %.4f)", est, lower, upper),
    pct_per_1C    = sprintf("%.2f%% (%.2f, %.2f)", pct_change, pct_lower, pct_upper)
  ) %>%
  dplyr::select(
    sex_label,
    baseline_rate,
    effect_per_1C,
    pct_per_1C
  )

print(summary_tbl)

