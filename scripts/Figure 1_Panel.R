source("scripts/functions.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(lfe)
  library(readr)
})

dir.create("outputs/raw_figures", recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# Load data
# ----------------------------
data_rok <- read.csv("inputs/SuicideData_ROK.csv")
data_rok$prec50 <- data_rok$prec / 50

# ----------------------------
# Main models
# ----------------------------
mod_temp <- felm(
  rate_adj ~ tmean + prec | fipsmo + stateyear | 0 | fips,
  data = data_rok,
  weights = data_rok$popw
)

mod_prec <- felm(
  rate_adj ~ tmean + prec50 | fipsmo + stateyear | 0 | fips,
  data = data_rok,
  weights = data_rok$popw
)

# ----------------------------
# Bootstrap results
# ----------------------------
boot_temp <- read.csv("inputs/bootstrap_runs/BootstrapMainModel_rok.csv")
boot_temp <- boot_temp[!is.na(boot_temp$est), ]

boot_prec <- read.csv("inputs/bootstrap_runs/BootstrapPrec50Model_rok.csv")
boot_prec <- boot_prec[!is.na(boot_prec$est), ]

# ----------------------------
# Panel (a): Temperature (straight line)
# ----------------------------
xx_temp  <- seq(-10, 30, by = 1)
ref_temp <- 12

beta_temp <- coef(mod_temp)["tmean"]
yy_temp   <- beta_temp * (xx_temp - ref_temp)

est_temp <- outer(boot_temp$est, (xx_temp - ref_temp))
ci_temp  <- apply(
  est_temp, 2,
  function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE)
)

ht_temp <- hist(data_rok$tmean, breaks = 40, plot = FALSE)

# ----------------------------
# Panel (b): Precipitation (straight line)
# ----------------------------
xx_prec_mm  <- seq(0, 250, by = 5)
xx_prec50   <- xx_prec_mm / 50
ref_prec_mm <- 100
ref_prec50  <- ref_prec_mm / 50

beta_prec <- coef(mod_prec)["prec50"]
yy_prec   <- beta_prec * (xx_prec50 - ref_prec50)

est_prec <- outer(boot_prec$est, (xx_prec50 - ref_prec50))
ci_prec  <- apply(
  est_prec, 2,
  function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE)
)

prec_trim <- data_rok$prec[data_rok$prec <= 250]
ht_prec <- hist(
  prec_trim,
  breaks = seq(0, 250, by = 10),
  plot = FALSE
)

# ----------------------------
# Helper: per-panel y scale
# ----------------------------
plot_with_hist <- function(
    x, y, ci, hist_obj,
    xlab, ylab, xlim, ref_x, panel_label
) {
  ymin <- min(ci, y, na.rm = TRUE)
  ymax <- max(ci, y, na.rm = TRUE)
  yr <- ymax - ymin
  if (yr < 0.03) yr <- 0.03
  
  ylim_use <- c(ymin - 0.22 * yr, ymax + 0.08 * yr)
  
  plot(
    x, y,
    type = "n",
    xlim = xlim,
    ylim = ylim_use,
    xlab = xlab,
    ylab = ylab,
    las = 1,
    cex.lab = 0.95,
    cex.axis = 0.9
  )
  
  polygon(
    c(x, rev(x)),
    c(ci[1, ], rev(ci[2, ])),
    col = "lightblue",
    border = NA
  )
  
  lines(x, y, lwd = 2)
  abline(h = 0, lty = 2, col = "black", lwd = 0.8)
  abline(v = ref_x, lty = 2, col = "black", lwd = 0.8)
  
  hist_base   <- ylim_use[1]
  hist_height <- 0.16 * diff(ylim_use)
  
  rect(
    xleft   = hist_obj$breaks[-length(hist_obj$breaks)],
    ybottom = hist_base,
    xright  = hist_obj$breaks[-1],
    ytop    = hist_base + (hist_obj$counts / max(hist_obj$counts)) * hist_height,
    col     = "lightblue",
    border  = "white"
  )
  
  mtext(panel_label, side = 3, adj = 0, line = 0.4, font = 2)
}

# ----------------------------
# Save figure
# ----------------------------
pdf("outputs/raw_figures/Figure_Temp_Prec_BothLinear.pdf", width = 12.5, height = 6.5)

par(
  mfrow = c(1, 2),
  mar = c(4.8, 5.2, 2.4, 1.2),
  mgp = c(2.8, 0.8, 0)
)

plot_with_hist(
  x = xx_temp,
  y = yy_temp,
  ci = ci_temp,
  hist_obj = ht_temp,
  xlab = "Monthly mean temperature (°C)",
  ylab = "Change in monthly suicide rate\n(relative to 12°C, per 100,000)",
  xlim = c(-10, 30),
  ref_x = ref_temp,
  panel_label = "(a) Temperature"
)

plot_with_hist(
  x = xx_prec_mm,
  y = yy_prec,
  ci = ci_prec,
  hist_obj = ht_prec,
  xlab = "Monthly precipitation (mm)",
  ylab = "Change in monthly suicide rate\n(relative to 100 mm, per 100,000)",
  xlim = c(0, 250),
  ref_x = ref_prec_mm,
  panel_label = "(b) Precipitation"
)

dev.off()