# functions_mod.R
# Backward-compatible helper file plus new utilities for heterogeneity models.
# Place this file in scripts/functions_mod.R, or source it from the working directory.

suppressPackageStartupMessages({
  library(dplyr)
  library(lfe)
})

# -----------------------------------------------------------------------------
# Original utility functions retained for backward compatibility
# -----------------------------------------------------------------------------

shift <- function(x, shift_by, y) {
  ana <- rep(NA, abs(shift_by))
  if (shift_by > 0) {
    out <- c(x[(1 + shift_by):length(x)], ana)
    ck <- c(y[(1 + shift_by):length(y)], ana)
    out[y != ck] <- NA
  } else if (shift_by < 0) {
    out <- c(ana, x[1:(length(x) + shift_by)])
    ck <- c(ana, y[1:(length(y) + shift_by)])
    out[y != ck] <- NA
  } else {
    out <- x
  }
  return(out)
}

"%&%" <- function(x, y) paste(x, y, sep = "")

fitmod <- function(data = dta, y = "rate_adj", vb = "tmean+prec",
                   fe = "fipsmo + stateyear", cl = "fips",
                   inter = "no", intvar = "", name = "") {
  if (inter == "yes") {
    data$above <- data[, intvar] > median(unlist(data[, intvar]), na.rm = TRUE)
    data$tmeanabove <- data$tmean * data$above
    data$precabove <- data$prec * data$above
    vb <- "tmean + tmeanabove + prec + precabove"
  }
  fmla <- as.formula(paste0(y, "~", vb, " | ", fe, " | 0 | ", cl))
  mod <- felm(fmla, data = data, weights = data$popw)
  coef_tbl <- summary(mod)$coefficients
  if (inter == "no") {
    mean_y <- weighted.mean(unlist(data[, y]), data$popw, na.rm = TRUE)
    out <- data.frame(
      est = coef_tbl["tmean", 1], se = coef_tbl["tmean", 2],
      mean = mean_y, var = name, group = "none", N = mod$N
    )
  } else {
    mns <- data %>%
      group_by(above) %>%
      summarise(mean = weighted.mean(rate_adj, popw, na.rm = TRUE), .groups = "drop")
    # Avoid importing multcomp here. Compute the linear combination directly.
    b <- coef(mod)
    V <- vcov(mod)
    est_high <- unname(b["tmean"] + b["tmeanabove"])
    se_high <- sqrt(V["tmean", "tmean"] + V["tmeanabove", "tmeanabove"] +
                      2 * V["tmean", "tmeanabove"])
    r1 <- data.frame(
      est = coef_tbl["tmean", 1], se = coef_tbl["tmean", 2],
      mean = mns$mean[mns$above == FALSE], var = name, group = "below", N = mod$N / 2
    )
    r2 <- data.frame(
      est = est_high, se = se_high,
      mean = mns$mean[mns$above == TRUE], var = name, group = "above", N = mod$N / 2
    )
    out <- rbind(r1, r2)
  }
  out
}

runreg <- function(lag = 0, lead = 0, y = "rate_adj", temp = "tmean", prec = "prec",
                   data = dta, weights = dta$popw, fe = "fipsmo + stateyear", cl = "fips") {
  if (lead == 0) {
    ll <- c(paste0(temp, "_lag", 0:lag), paste0(prec, "_lag", 0:lag))
  } else {
    ll <- c(
      paste0(temp, "_lead", lead:1), paste0(temp, "_lag", 0:lag),
      paste0(prec, "_lead", lead:1), paste0(prec, "_lag", 0:lag)
    )
  }
  lgs <- paste(ll, collapse = "+")
  fmla <- as.formula(paste0(y, " ~", lgs, " | ", fe, " | 0 | ", cl))
  mod <- felm(fmla, data = data, weights = weights)
  coef_tbl <- summary(mod)$coefficients
  mn <- weighted.mean(unlist(data[, y]), weights, na.rm = TRUE)
  vars <- 1:(lag + lead + 1)
  vars1 <- (1 + lead):(lag + lead + 1)
  r1 <- data.frame(
    est = coef_tbl[vars, 1], se = coef_tbl[vars, 2], mean = mn,
    var = "leadlag", group = c("t+1", "t", "t-1")[seq_along(vars)], N = mod$N
  )
  r2 <- data.frame(
    est = sum(mod$coefficients[vars1]),
    se = sqrt(sum(vcov(mod)[vars1, vars1])),
    mean = mn, var = "leadlag", group = "combined", N = mod$N
  )
  rbind(r1, r2)
}

run_reg_bstrap <- function(data, clustername) {
  sub_data <- getdata(data, cluster(data = data, clustername = clustername,
                                    size = n_clust, method = "srswr"))
  output <- sum(summary(felm(fmla, weights = sub_data$popw, data = sub_data))$coefficients[c("tmean", "tmean_lag1"), "Estimate"])
  return(output)
}

# -----------------------------------------------------------------------------
# New helper functions for revised heterogeneity models
# -----------------------------------------------------------------------------

weighted_center <- function(x, w = NULL) {
  if (is.null(w)) {
    center <- mean(x, na.rm = TRUE)
  } else {
    center <- weighted.mean(x, w, na.rm = TRUE)
  }
  x - center
}

get_term <- function(coefs, candidates) {
  hit <- candidates[candidates %in% names(coefs)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

lincom_effect <- function(model, base_term, interaction_term = NULL, x = 0, scale = 1) {
  b <- coef(model)
  V <- vcov(model)
  if (!base_term %in% names(b)) {
    return(data.frame(est = NA_real_, se = NA_real_, lower = NA_real_, upper = NA_real_))
  }
  est <- unname(b[base_term])
  var <- V[base_term, base_term]
  if (!is.null(interaction_term) && !is.na(interaction_term) && interaction_term %in% names(b)) {
    est <- est + x * unname(b[interaction_term])
    var <- var + x^2 * V[interaction_term, interaction_term] + 2 * x * V[base_term, interaction_term]
  }
  est <- est * scale
  se <- sqrt(var) * abs(scale)
  data.frame(est = est, se = se, lower = est - 1.96 * se, upper = est + 1.96 * se)
}

coef_extract <- function(model, term, label = term, scale = 1) {
  b <- coef(model)
  V <- vcov(model)
  if (!term %in% names(b)) {
    return(data.frame(term = label, est = NA_real_, se = NA_real_, t = NA_real_, p = NA_real_))
  }
  se <- sqrt(V[term, term])
  est <- unname(b[term])
  tval <- est / se
  pval <- 2 * pnorm(abs(tval), lower.tail = FALSE)
  data.frame(term = label, est = est * scale, se = se * abs(scale), t = tval, p = pval)
}
