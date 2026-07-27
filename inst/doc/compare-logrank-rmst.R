## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align = "center",
  fig.width = 7,
  fig.height = 4.2
)
library(FastSurvival)
has_gsDesign <- requireNamespace("gsDesign", quietly = TRUE)
knitr::opts_chunk$set(eval = has_gsDesign)
nsim <- 10000
set.seed(1)

## ----gsdesign-note, eval = !has_gsDesign, echo = FALSE, results = "asis"------
# cat("This vignette computes the sample size with the gsDesign package.",
#     "Install gsDesign to evaluate the code and reproduce the figures and the",
#     "power table below.")

## ----design-------------------------------------------------------------------
m0             <- 12
lam0           <- log(2) / m0
hr_ph          <- 0.75
hr_late        <- 0.60
hr_early_cross <- 1.40
delay          <- 6
alpha          <- 0.025
power_target   <- 0.90
accrual        <- 12
minfup         <- 36
study_dur      <- accrual + minfup
dropout_annual <- 0.05
eta            <- -log(1 - dropout_annual) / 12
tau            <- 30

ns <- gsDesign::nSurv(
  lambdaC = lam0, hr = hr_ph, eta = eta,
  T = study_dur, minfup = minfup,
  alpha = alpha, beta = 1 - power_target, sided = 1, ratio = 1
)
n_per   <- ceiling(ns$n / 2)
n       <- c(n_per, n_per)
n_total <- sum(n)
a_rate  <- n_total / accrual

## ----design-print-------------------------------------------------------------
data.frame(
  "Total n"       = n_total,
  "Per group"     = n_per,
  "Target events" = round(ns$d),
  check.names = FALSE
)

## ----scenarios, fig.width = 9, fig.height = 3.6-------------------------------
scn <- gen_scenario_fast(
  scenarios = list(
    "PH" = list(
      e.hazard = list(lam0, hr_ph * lam0)
    ),
    "Delayed effect" = list(
      e.hazard = list(lam0, c(lam0, hr_late * lam0)),
      e.time   = c(0, delay, Inf)
    ),
    "Crossing" = list(
      e.hazard = list(lam0, c(hr_early_cross * lam0, hr_late * lam0)),
      e.time   = c(0, delay, Inf)
    )
  ),
  shared = list(n = n, a.time = c(0, accrual), a.rate = a_rate)
)

plot(scn, tmax = study_dur, mfrow = c(1, 3))

## ----simulate-----------------------------------------------------------------
scenarios <- scn$scenarios
seeds     <- c(101, 102, 103)

power_tab <- data.frame(
  Scenario = character(0), LogRank = numeric(0), RMST = numeric(0),
  stringsAsFactors = FALSE
)
examples <- vector("list", length(scenarios))
names(examples) <- names(scenarios)

for (i in seq_along(scenarios)) {
  s   <- scenarios[[i]]
  dat <- do.call(
    simdata_fast,
    c(s$args, list(nsim = nsim, d.hazard = eta, seed = seeds[i]))
  )
  res <- analysis_fast(
    dat, control = 1,
    time.looks = study_dur,
    stat = c("logrank", "rmst"),
    tau = tau, side = 1
  )
  s_lr   <- simsummary_fast(res, p.col = "logrank.p", alpha = alpha)
  s_rmst <- simsummary_fast(res, p.col = "rmst.p",    alpha = alpha)
  power_tab <- rbind(power_tab, data.frame(
    Scenario = s$label,
    LogRank  = s_lr[s_lr$look   == "overall", "cum.reject"],
    RMST     = s_rmst[s_rmst$look == "overall", "cum.reject"],
    stringsAsFactors = FALSE
  ))
  examples[[i]] <- dat[dat$sim == 1L, c("tte", "event", "group")]
}

## ----km-plot, fig.width = 9, fig.height = 4.6---------------------------------
ex  <- examples[["Crossing"]]
fit <- kmcurve_fast(ex$tte, ex$event, ex$group, control = 1)
plot(fit, hr = TRUE, rmst = TRUE, tau = tau, bw = 3)

## ----power-table--------------------------------------------------------------
knitr::kable(
  power_tab, digits = 3,
  col.names = c("Scenario", "Log-rank", "RMST"),
  caption = "Simulated power at the fixed analysis (one-sided 0.025)."
)

