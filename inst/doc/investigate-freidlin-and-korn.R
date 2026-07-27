## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 9,
                      fig.height = 4, dpi = 110)

## ----library------------------------------------------------------------------
library(FastSurvival)

## ----scenarios----------------------------------------------------------------
e_time <- c(0, 0.1, Inf)         # change point at 0.1 year
cut    <- 5                      # analysis at 5 years

scenarios <- list(
  equal = list(
    label = "Equal survival (weak null)",
    haz_c = c(0.25, 0.25), haz_e = c(0.25, 0.25)),
  strong = list(
    label = "Experimental uniformly worse (strong null)",
    haz_c = c(0.25, 0.25), haz_e = c(16, 0.76) * 0.25)
)

## ----scenario-figure, echo = FALSE, fig.width = 5, fig.height = 4-------------
pe_cumhaz <- function(t, haz, e_time) {
  ch <- numeric(length(t))
  for (j in seq_along(haz)) {
    ch <- ch + haz[j] * pmax(0, pmin(t, e_time[j + 1]) - e_time[j])
  }
  ch
}
pe_surv <- function(t, haz, e_time) exp(-pe_cumhaz(t, haz, e_time))

col_control <- "#1b9e77"
col_treat   <- "#d95f02"
scn <- scenarios$strong
tt  <- seq(0, cut, length.out = 601)

plot(tt, pe_surv(tt, scn$haz_c, e_time), type = "l", col = col_control, lwd = 2,
     ylim = c(0, 1), xlab = "Time (years)", ylab = "Survival probability",
     main = "Survival (strong null)")
lines(tt, pe_surv(tt, scn$haz_e, e_time), col = col_treat, lwd = 2)
legend("topright", bty = "n", cex = 0.85,
       legend = c("Control", "Experimental"),
       col = c(col_control, col_treat), lwd = 2)

## ----design-------------------------------------------------------------------
NSIM    <- 2000L     # raise (e.g. 5000) for tighter Monte Carlo error
N_TOTAL <- 2000L
ALPHA   <- 0.025
SIDE    <- 1L
T_STAR  <- 0.5       # modestly-weighted log-rank delay (years)
SEED    <- 20240601L

run_one_scenario <- function(scn) {
  dataset <- simdata_fast(
    nsim     = NSIM,
    n        = rep(N_TOTAL %/% 2L, 2L),
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
  res_mc   <- analysis_fast(dataset, control = 1, time.looks = cut,
                            stat = "maxcombo",
                            mc.rho = c(0, 0, 1, 1), mc.gamma = c(0, 1, 0, 1),
                            side = SIDE)

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
round(rates, 3)

## ----results-table, echo = FALSE----------------------------------------------
knitr::kable(round(rates, 3),
             caption = paste0("One-sided rejection rate in favor of the ",
                              "experimental arm (", NSIM, " simulated trials, ",
                              "n = ", N_TOTAL, "). Under the weak null every ",
                              "valid test is near 2.5%; under the strong null ",
                              "the experimental arm is uniformly worse, so a ",
                              "rate above 2.5% reflects a spurious finding."))

