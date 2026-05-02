suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(lfe)
  library(ggplot2)
  library(tibble)
})

# =========================================================
# 0) 경로 설정
# =========================================================
base_path <- "inputs"
out_dir   <- file.path(base_path, "outputs_age")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

f_base <- file.path(base_path, "SuicideData_ROK.csv")
f_age  <- file.path(base_path, "suicide_age.csv")
f_pop  <- file.path(base_path, "age_population.csv")

stopifnot(file.exists(f_base), file.exists(f_age), file.exists(f_pop))

# =========================================================
# 1) 데이터 불러오기
# =========================================================
base <- readr::read_csv(f_base, show_col_types = FALSE)
suD  <- readr::read_csv(f_age,  show_col_types = FALSE)
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

suD <- suD %>%
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
# 4) suicide_age long 변환
#    age_3 ~ age_19 사용
#    age code:
#      3=5-9, 4=10-14, ..., 15=65-69, ..., 19=85+
# =========================================================
age_vars <- paste0("age_", 3:19)

missing_age_vars <- setdiff(age_vars, names(suD))
if (length(missing_age_vars) > 0) {
  stop("suicide_age.csv에 다음 변수가 없습니다: ",
       paste(missing_age_vars, collapse = ", "))
}

su_long <- suD %>%
  dplyr::select(fips, year, month, dplyr::all_of(age_vars)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(age_vars),
    names_to = "age_var",
    values_to = "deaths"
  ) %>%
  dplyr::mutate(
    age_code = as.integer(stringr::str_remove(age_var, "age_"))
  )

# =========================================================
# 5) 연령코드 -> 연령대 라벨 / midpoint 생성
#    midpoint는 상호작용 회귀용
# =========================================================
age_map <- tibble::tibble(
  age_code  = 3:19,
  age_label = c(
    "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39",
    "40-44", "45-49", "50-54", "55-59", "60-64", "65-69",
    "70-74", "75-79", "80-84", "85+"
  ),
  age_mid = c(7, 12, 17, 22, 27, 32, 37, 42, 47, 52, 57, 62, 67, 72, 77, 82, 87)
)

su_long <- su_long %>%
  dplyr::left_join(age_map, by = "age_code")

# =========================================================
# 6) age_population 처리
#    sex==3 (전체)만 사용
#    age code 3:19 사용
#    pop_1997 ~ pop_2024 -> long
# =========================================================
pop_long <- popD %>%
  dplyr::filter(sex == 3, age %in% 3:19) %>%
  tidyr::pivot_longer(
    cols = dplyr::starts_with("pop_"),
    names_to = "year_var",
    values_to = "population"
  ) %>%
  dplyr::mutate(
    year = as.integer(stringr::str_remove(year_var, "pop_")),
    age_code = age
  ) %>%
  dplyr::select(fips, year, age_code, population) %>%
  tidyr::crossing(month = 1:12) %>%
  dplyr::mutate(month = as.integer(month)) %>%
  dplyr::left_join(age_map, by = "age_code")

# =========================================================
# 7) 병합 및 자살률 계산
# =========================================================
dat_age <- su_long %>%
  dplyr::left_join(
    pop_long,
    by = c("fips", "year", "month", "age_code", "age_label", "age_mid")
  ) %>%
  dplyr::left_join(
    base2,
    by = c("fips", "year", "month")
  ) %>%
  dplyr::mutate(
    deaths     = as.numeric(deaths),
    population = as.numeric(population),
    rate_adj   = dplyr::if_else(
      !is.na(population) & population > 0,
      deaths / population * 100000,
      NA_real_
    )
  ) %>%
  dplyr::filter(!is.na(rate_adj), !is.na(tmean), !is.na(prec), !is.na(age_mid))

cat("Final analytic sample n =", nrow(dat_age), "\n")
cat("Year range:", min(dat_age$year, na.rm = TRUE), "~", max(dat_age$year, na.rm = TRUE), "\n")

# =========================================================
# 8) 연령 변수 중심화
#    5년 단위 해석을 쉽게 하기 위해 (age_mid - mean)/5
# =========================================================
age_mid_mean <- mean(dat_age$age_mid, na.rm = TRUE)

dat_age <- dat_age %>%
  dplyr::mutate(
    age5_c = (age_mid - age_mid_mean) / 5
  )

cat("Mean age midpoint =", round(age_mid_mean, 3), "\n")

# =========================================================
# 9) baseline 자살률 (연령대별)
# =========================================================
baseline_tbl <- dat_age %>%
  dplyr::group_by(age_code, age_label, age_mid) %>%
  dplyr::summarise(
    total_deaths   = sum(deaths, na.rm = TRUE),
    total_pop      = sum(population, na.rm = TRUE),
    baseline_rate  = total_deaths / total_pop * 100000,
    mean_rate_cell = mean(rate_adj, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(age_code)

print(baseline_tbl, n = Inf)

# =========================================================
# 10) 연령 상호작용 회귀
#     age5_c 1단위 = 5세 증가
# =========================================================
mod_age <- lfe::felm(
  rate_adj ~ tmean + age5_c + tmean:age5_c + prec |
    fipsmo + stateyear | 0 | fips,
  data = dat_age,
  weights = dat_age$population
)

mod_sum <- summary(mod_age)
print(mod_sum)

coef_df <- as.data.frame(mod_sum$coefficients)
coef_df$term <- rownames(coef_df)
rownames(coef_df) <- NULL

# =========================================================
# 11) 연령별 1℃ 효과 계산
#     effect(age) = b1 + b3 * age5_c
# =========================================================
b <- coef(mod_age)
V <- vcov(mod_age)

int_name <- if ("tmean:age5_c" %in% names(b)) "tmean:age5_c" else "age5_c:tmean"

b1 <- unname(b["tmean"])
b3 <- unname(b[int_name])

age_effect_tbl <- age_map %>%
  dplyr::mutate(
    age5_c = (age_mid - age_mid_mean) / 5,
    est = b1 + b3 * age5_c,
    var = V["tmean", "tmean"] +
      (age5_c^2) * V[int_name, int_name] +
      2 * age5_c * V["tmean", int_name],
    se = sqrt(var),
    lower = est - 1.96 * se,
    upper = est + 1.96 * se
  ) %>%
  dplyr::left_join(
    baseline_tbl %>% dplyr::select(age_code, baseline_rate),
    by = "age_code"
  ) %>%
  dplyr::mutate(
    pct_change = est / baseline_rate * 100,
    pct_lower  = lower / baseline_rate * 100,
    pct_upper  = upper / baseline_rate * 100
  ) %>%
  dplyr::arrange(age_code)

print(age_effect_tbl, n = Inf)

# =========================================================
# 12) 대표 연령군 비교표 (예: 20-39 vs 65+)
# =========================================================
summary_groups <- age_effect_tbl %>%
  dplyr::mutate(
    age_group = dplyr::case_when(
      age_mid >= 22 & age_mid <= 37 ~ "20-39",
      age_mid >= 67                ~ "65+",
      TRUE                         ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(age_group)) %>%
  dplyr::group_by(age_group) %>%
  dplyr::summarise(
    baseline_rate = mean(baseline_rate, na.rm = TRUE),
    est = mean(est, na.rm = TRUE),
    lower = mean(lower, na.rm = TRUE),
    upper = mean(upper, na.rm = TRUE),
    pct_change = mean(pct_change, na.rm = TRUE),
    pct_lower = mean(pct_lower, na.rm = TRUE),
    pct_upper = mean(pct_upper, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    effect_per_1C = sprintf("%.4f (%.4f, %.4f)", est, lower, upper),
    pct_per_1C    = sprintf("%.2f%% (%.2f, %.2f)", pct_change, pct_lower, pct_upper)
  )

print(summary_groups)

