suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(lfe)
  library(ggplot2)
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
    stateyear = floor(fips / 1000) * 10000 + year
  )

# -----------------------------
# 2) grdp wide -> long
# -----------------------------
grdp <- fread("inputs/grdp.csv") %>% as.data.frame()
names(grdp)[1] <- "fips"

grdp_long <- tidyr::pivot_longer(
  grdp,
  cols = -fips,
  names_to = "year",
  values_to = "grdp"
)

grdp_long <- grdp_long %>%
  mutate(
    fips = as.numeric(fips),
    year = as.integer(year),
    grdp = as.numeric(grdp),
    ln_grdp = log(grdp)
  ) %>%
  select(fips, year, grdp, ln_grdp)

# -----------------------------
# 3) merge
# -----------------------------
dat2 <- dat %>%
  left_join(grdp_long, by = c("fips", "year"))

cat("총 관측치 수:", nrow(dat2), "\n")
cat("ln_grdp 결측 수:", sum(is.na(dat2$ln_grdp)), "\n")

# -----------------------------
# 4) regression
# -----------------------------
mod_income <- felm(
  rate_adj ~ tmean + ln_grdp + tmean:ln_grdp + prec |
    fipsmo + stateyear |
    0 |
    fips,
  data = dat2,
  weights = dat2$popw
)

summary(mod_income)

# -----------------------------
# 5) marginal effect 계산
# -----------------------------
b <- coef(mod_income)
V <- vcov(mod_income)

int_name <- if ("tmean:ln_grdp" %in% names(b)) {
  "tmean:ln_grdp"
} else if ("ln_grdp:tmean" %in% names(b)) {
  "ln_grdp:tmean"
} else {
  stop("상호작용항 이름을 찾을 수 없습니다.")
}

b1 <- b["tmean"]
b3 <- b[int_name]

grid <- seq(
  quantile(dat2$ln_grdp, 0.05, na.rm = TRUE),
  quantile(dat2$ln_grdp, 0.95, na.rm = TRUE),
  length.out = 100
)

me <- b1 + b3 * grid

se <- sapply(grid, function(g) {
  sqrt(
    V["tmean", "tmean"] +
      g^2 * V[int_name, int_name] +
      2 * g * V["tmean", int_name]
  )
})

df_me <- data.frame(
  ln_grdp = grid,
  me = me,
  lo = me - 1.96 * se,
  hi = me + 1.96 * se
)

# -----------------------------
# 6) p10 / p90 효과 및 차이 계산
# -----------------------------
p10 <- quantile(dat2$ln_grdp, 0.10, na.rm = TRUE)
p90 <- quantile(dat2$ln_grdp, 0.90, na.rm = TRUE)

me_p10 <- b1 + b3 * p10
me_p90 <- b1 + b3 * p90

delta <- me_p90 - me_p10
se_delta <- abs(p90 - p10) * sqrt(V[int_name, int_name])
ci_low <- delta - 1.96 * se_delta
ci_high <- delta + 1.96 * se_delta

# -----------------------------
# 7) plot
# -----------------------------
label_df <- data.frame(
  x = c(p10, p90),
  y = c(me_p10, me_p90),
  lab = c("Bottom 10%", "Top 10%")
)

x_annot <- p90 - 0.15
y_annot <- max(df_me$hi, na.rm = TRUE) * 0.92

ggplot(df_me, aes(x = ln_grdp, y = me)) +
  geom_ribbon(aes(ymin = lo, ymax = hi),
              fill = "lightblue", alpha = 0.35) +
  geom_line(linewidth = 1.1, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = c(p10, p90), linetype = "dotted", color = "gray35") +
  geom_point(aes(x = p10, y = me_p10), size = 3, color = "black") +
  geom_point(aes(x = p90, y = me_p90), size = 3, color = "black") +
  geom_text(
    data = label_df,
    aes(x = x, y = y, label = lab),
    nudge_y = 0.01,
    size = 4.2,
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = x_annot,
    y = y_annot,
    label = paste0(
      "\u0394 (Top 10% - Bottom 10%) = ",
      round(delta, 3),
      "\n95% CI [",
      round(ci_low, 3), ", ", round(ci_high, 3), "]"
    ),
    hjust = 1,
    size = 4.2,
    fontface = "bold"
  ) +
  labs(
    x = "Log GRDP",
    y = "Marginal effect of temperature on suicide rate",
    title = "Marginal Effect of Temperature by Income Level",
    subtitle = "Shaded area indicates 95% confidence interval"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(size = 12),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )