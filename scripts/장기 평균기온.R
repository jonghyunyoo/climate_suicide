library(dplyr)
library(lfe)

dat <- read.csv("inputs/SuicideData_ROK.csv")
baseline_temp <- dat %>%
  group_by(fips) %>%
  summarise(
    avg_temp = mean(tmean, na.rm = TRUE),
    .groups = "drop"
  )
dat <- dat %>%
  left_join(baseline_temp, by = "fips")

mod_cont <- felm(
  rate_adj ~ tmean + avg_temp + tmean:avg_temp + prec |
    fipsmo + stateyear |
    0 |
    fips,
  data = dat,
  weights = dat$popw
)

summary(mod_cont)

dat <- dat %>%
  mutate(
    warm_region = ifelse(avg_temp >= 12, 1, 0)
  )

mod_bin <- felm(
  rate_adj ~ tmean + warm_region + tmean:warm_region + prec |
    fipsmo + stateyear |
    0 |
    fips,
  data = dat,
  weights = dat$popw
)

summary(mod_bin)