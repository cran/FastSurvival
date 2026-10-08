## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>"
)
has_gsd <- requireNamespace("gsDesign", quietly = TRUE)

## ----load---------------------------------------------------------------------
library(FastSurvival)

## ----simulate-----------------------------------------------------------------
nsim <- 10000
eta  <- -log(1 - 0.05) / 12
df <- simdata_fast(
  nsim     = nsim,
  n        = c(241, 241),
  a.time   = c(0, 23),
  a.prop   = 1,
  e.median = list(12.9, 9.0),
  d.hazard = list(eta, eta),
  seed     = 1
)

## ----analyze------------------------------------------------------------------
res <- analysis_fast(
  df, control = 2,
  event.looks = c(252, 336),
  stat = c("logrank", "coxph"), side = 2
)

## ----boundaries, eval = has_gsd-----------------------------------------------
gsd <- gsDesign::gsDesign(
  k         = 2,
  timing    = c(252, 336) / 336,
  alpha     = 0.025,
  beta      = 0.1,
  sfu       = gsDesign::sfLDOF,
  test.type = 1
)

spend_alpha <- 2 * pnorm(gsd$upper$bound, lower.tail = FALSE)
spend_alpha

## ----summary, eval = has_gsd--------------------------------------------------
oc <- simsummary_fast(res, p.col = "logrank.p", alpha = spend_alpha)
oc

## ----closed-form, eval = has_gsd----------------------------------------------
looks <- c(252, 336)
theta <- log(12.9 / 9.0) / 2
gp <- gsDesign::gsProbability(k = 2, theta = theta, n.I = looks,
                              a = rep(-20, 2), b = gsd$upper$bound)
p_cross <- gp$upper$prob[, 1]

# Expected number of deaths by calendar time t: uniform accrual of 482
# subjects over 23 months, exponential survival and dropout in each group.
lam <- log(2) / c(12.9, 9.0)
exp_events <- function(t, n_tot = 482, r_acc = 23) {
  a_end <- min(t, r_acc)
  sum(vapply(lam, function(l) {
    cc <- l + eta
    (n_tot / 2) / r_acc * l / cc *
      (a_end - (exp(-cc * (t - a_end)) - exp(-cc * t)) / cc)
  }, numeric(1)))
}
t_look <- vapply(looks, function(d) {
  stats::uniroot(function(t) exp_events(t) - d, c(1, 200))$root
}, numeric(1))

sim_val <- function(lk, col) as.numeric(oc[oc$look == lk, col])
z1 <- -res$logrank.z[res$look == 1]
cmp <- data.frame(
  quantity = c("Mean log-rank Z at the interim",
               "Probability of crossing at the interim",
               "Overall power",
               "Expected number of events at stop",
               "Expected analysis time at stop (months)"),
  simulation = c(mean(z1),
                 sim_val("1", "prob.stop.efficacy"),
                 sim_val("overall", "cum.reject"),
                 sim_val("overall", "n.event.mean"),
                 sim_val("overall", "cutoff.mean")),
  closed_form = c(theta * sqrt(looks[1]),
                  p_cross[1],
                  sum(p_cross),
                  sum(looks * c(p_cross[1], 1 - p_cross[1])),
                  sum(t_look * c(p_cross[1], 1 - p_cross[1])))
)
knitr::kable(cmp, digits = 3,
             col.names = c("Quantity",
                           paste0("FastSurvival (", nsim, " trials)"),
                           "Closed form"))

## ----nph, eval = has_gsd------------------------------------------------------
nsim_nph <- 2000
df_delay <- simdata_fast(
  nsim     = nsim_nph,
  n        = c(241, 241),
  a.time   = c(0, 23),
  a.prop   = 1,
  e.hazard = list(c(0.077, 0.045), c(0.077, 0.077)),
  e.time   = c(0, 3, Inf),
  d.hazard = list(eta, eta),
  seed     = 1
)

set.seed(1)
res_delay <- analysis_fast(
  df_delay, control = 2,
  event.looks = looks,
  stat = c("logrank", "maxcombo"), side = 2,
  mc.alpha = spend_alpha
)

oc_lr <- simsummary_fast(res_delay, p.col = "logrank.p",  alpha = spend_alpha)
oc_mc <- simsummary_fast(res_delay, p.col = "maxcombo.p", alpha = spend_alpha)
data.frame(
  test  = c("Log-rank", "Max-combo"),
  power = round(c(oc_lr[oc_lr$look == "overall", "cum.reject"],
                  oc_mc[oc_mc$look == "overall", "cum.reject"]), 3)
)

## ----nph-share, eval = has_gsd, echo = FALSE, results = "asis"----------------
cat("Of the ", format(2 * nsim_nph, big.mark = ","),
    " max-combo p-values (two looks per simulated trial), ",
    round(100 * mean(res_delay$maxcombo.p.exact, na.rm = TRUE)),
    "% required the multivariate normal integral.\n", sep = "")

