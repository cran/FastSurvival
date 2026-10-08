## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 9,
                      fig.height = 4, dpi = 110)

## ----library------------------------------------------------------------------
library(FastSurvival)

## ----scenarios----------------------------------------------------------------
e_time <- c(0, 0.1, Inf)         # change point at 0.1 year
cut    <- 5                      # analysis at 5 years
N_HALF <- 1000L                  # subjects per arm

scenarios <- list(
  equal = list(
    label = "Equal survival (weak null)",
    haz_c = c(0.25, 0.25), haz_e = c(0.25, 0.25)),
  strong = list(
    label = "Experimental uniformly worse (strong null)",
    haz_c = c(0.25, 0.25), haz_e = c(16, 0.76) * 0.25)
)

## ----scenario-figure, fig.width = 5, fig.height = 4---------------------------
scn_strong <- gen_scenario_fast(
  list("Strong null" = list(
    e.hazard = list(scenarios$strong$haz_c, scenarios$strong$haz_e),
    e.time   = e_time)),
  shared = list(n = c(N_HALF, N_HALF), a.time = c(0, 1e-4), a.prop = 1)
)
plot(scn_strong, tmax = cut, hr_max = 2, xlab = "Time (years)")

## ----design-------------------------------------------------------------------
NSIM    <- 2000L     # raise (e.g. 5000) for tighter Monte Carlo error
N_TOTAL <- 2L * N_HALF
ALPHA   <- 0.025
SIDE    <- 1L
T_STAR  <- 0.5       # modestly-weighted log-rank delay (years)
SEED    <- 20240601L

run_one_scenario <- function(scn) {
  dataset <- simdata_fast(
    nsim     = NSIM,
    n        = c(N_HALF, N_HALF),
    a.time   = c(0, 1e-4),       # near-instantaneous enrollment
    a.prop   = 1,
    e.hazard = list(scn$haz_c, scn$haz_e),
    e.time   = e_time,
    seed     = SEED
  )

  res_lr   <- analysis_fast(dataset, control = 1, time.looks = cut,
                            stat = "logrank", side = SIDE)
  res_fh01 <- analysis_fast(dataset, control = 1, time.looks = cut,
                            stat = "logrank", weight = "fh",
                            rho = 0, gamma = 1, side = SIDE)
  res_mb   <- analysis_fast(dataset, control = 1, time.looks = cut,
                            stat = "logrank", weight = "mwlrt",
                            t_star = T_STAR, side = SIDE)
  set.seed(SEED)
  res_mc   <- analysis_fast(dataset, control = 1, time.looks = cut,
                            stat = "maxcombo",
                            mc.rho = c(0, 0, 1, 1), mc.gamma = c(0, 1, 0, 1),
                            side = SIDE, mc.alpha = ALPHA)

  reject_in_favor <- function(res, pcol) {
    ok <- res$reached
    mean(res[[pcol]][ok] <= ALPHA, na.rm = TRUE)
  }
  c("Log-rank"          = reject_in_favor(res_lr,   "logrank.p"),
    "FH(0,1)"           = reject_in_favor(res_fh01, "logrank.p"),
    "Modestly weighted" = reject_in_favor(res_mb,   "logrank.p"),
    "MaxCombo"          = reject_in_favor(res_mc,   "maxcombo.p"))
}

## ----run-all------------------------------------------------------------------
rates <- vapply(scenarios, run_one_scenario, numeric(4L))
colnames(rates) <- vapply(scenarios, function(s) s$label, character(1L))

## ----results-table, echo = FALSE----------------------------------------------
knitr::kable(round(rates, 3),
             caption = paste0("One-sided rejection rate in favor of the ",
                              "experimental arm (", NSIM, " simulated trials, ",
                              "n = ", N_TOTAL, "). Under the weak null every ",
                              "valid test is near 2.5%; under the strong null ",
                              "the experimental arm is uniformly worse, so a ",
                              "rate above 2.5% reflects a spurious finding."))

