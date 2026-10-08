## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6.5,
  fig.height = 4.2,
  fig.align = "center"
)
library(FastSurvival)

## ----design-------------------------------------------------------------------
nsim <- 2000
seed <- 20260711

# Three arms: group 1 = control, groups 2 and 3 = experimental.
# Per-arm sample size (balanced 1:1:1 randomization).
n_arm <- c(150, 150, 150)

# Median survival by arm, in months. The two experimental arms have an
# increasing benefit over control.
mst <- c(control = 12, dose_low = 16, dose_high = 20)

# Piecewise-uniform accrual: constant rate over 12 months enrolling all 450.
a_time <- c(0, 12)
a_rate <- sum(n_arm) / 12

# Single calendar analysis time, in months.
analysis_month <- 30

# One-sided family-wise significance level and the number of experimental arms.
alpha_total <- 0.025
n_contrast  <- length(n_arm) - 1L

# Implied exponential hazard ratios of each experimental arm versus control.
data.frame(
  arm = names(mst)[-1],
  median = mst[-1],
  HR_vs_control = round(mst["control"] / mst[-1], 3),
  row.names = NULL
)

## ----generate-----------------------------------------------------------------
df <- simdata_fast(
  nsim     = nsim,
  n        = n_arm,
  a.time   = a_time,
  a.rate   = a_rate,
  e.median = list(mst["control"], mst["dose_low"], mst["dose_high"]),
  seed     = seed
)

# One row per subject per simulation; the group column labels the arm.
table(df$group[df$sim == 1])

## ----pairwise-----------------------------------------------------------------
control_arm <- 1L
exp_arms    <- setdiff(sort(unique(df$group)), control_arm)

# Matrix of one-sided log-rank statistics and p-values, one column per contrast.
Z <- matrix(NA_real_, nrow = nsim, ncol = n_contrast)
P <- matrix(NA_real_, nrow = nsim, ncol = n_contrast)
colnames(Z) <- colnames(P) <- paste0("arm", exp_arms, "_vs_ctrl")

for (k in seq_along(exp_arms)) {
  sub <- df[df$group %in% c(control_arm, exp_arms[k]), ]
  res <- analysis_fast(
    sub,
    control    = control_arm,
    time.looks = analysis_month,
    stat       = "logrank",
    side       = 1
  )
  Z[, k] <- res$logrank.z
  P[, k] <- res$logrank.p
}

head(round(P, 4))

## ----pairwise-fast------------------------------------------------------------
pw <- pairwise_fast(
  df,
  control    = control_arm,
  time.looks = analysis_month,
  stat       = "logrank",
  side       = 1,
  adjust     = "bonferroni"
)
head(pw[, c("arm", "sim", "cutoff", "logrank.z", "logrank.p", "p.adj")])

## ----pairwise-fast-check------------------------------------------------------
round(tapply(pw$logrank.p <= alpha_total, pw$arm, mean), 3)

## ----pairwise-fast-event------------------------------------------------------
pw_ev <- pairwise_fast(
  df,
  control     = control_arm,
  event.looks = 200,
  primary     = 3,
  stat        = "logrank",
  side        = 1,
  adjust      = "bonferroni"
)
head(pw_ev[, c("arm", "sim", "cutoff", "reached", "logrank.z", "p.adj")])

## ----oc-alt-------------------------------------------------------------------
thr_marginal   <- alpha_total
thr_bonferroni <- alpha_total / n_contrast

reject_marginal   <- P <= thr_marginal
reject_bonferroni <- P <= thr_bonferroni

power_table <- data.frame(
  contrast          = colnames(P),
  power_marginal    = round(colMeans(reject_marginal), 3),
  power_bonferroni  = round(colMeans(reject_bonferroni), 3),
  row.names = NULL
)
power_table

## ----oc-disjunctive-----------------------------------------------------------
any_reject   <- rowSums(reject_bonferroni) > 0
# The stronger experimental arm is the last column (largest median).
stronger_col <- n_contrast
select_stronger <- reject_bonferroni[, stronger_col]

data.frame(
  quantity = c("Disjunctive power (any arm beats control)",
               "Power to declare the stronger arm",
               "Expected number of arms declared"),
  value = c(round(mean(any_reject), 3),
            round(mean(select_stronger), 3),
            round(mean(rowSums(reject_bonferroni)), 3)),
  row.names = NULL
)

## ----correlation--------------------------------------------------------------
round(stats::cor(Z, use = "complete.obs"), 3)

## ----null-fwer----------------------------------------------------------------
df0 <- simdata_fast(
  nsim     = nsim,
  n        = n_arm,
  a.time   = a_time,
  a.rate   = a_rate,
  e.median = list(mst["control"], mst["control"], mst["control"]),
  seed     = seed + 1
)

pw0 <- pairwise_fast(df0, control = control_arm,
                     time.looks = analysis_month, stat = "logrank", side = 1)
# One column per contrast; the rows of each arm are ordered by simulation.
P0 <- sapply(exp_arms, function(k) pw0$logrank.p[pw0$arm == k])

pce_hat  <- mean(P0 <= alpha_total)
fwer_hat <- mean(rowSums(P0 <= thr_bonferroni) > 0)
# Normal approximation with the shared-control correlation of one half.
fwer_nrm <- 1 - mvtnorm::pmvnorm(
  upper = rep(stats::qnorm(1 - thr_bonferroni), n_contrast),
  corr = matrix(0.5, n_contrast, n_contrast) + diag(0.5, n_contrast),
  algorithm = mvtnorm::TVPACK()
)[1]

data.frame(
  quantity = c("Per-contrast type I error (nominal 0.025)",
               "Family-wise error, Bonferroni",
               "Family-wise error, Bonferroni (normal approximation)"),
  estimate = round(c(pce_hat, fwer_hat, fwer_nrm), 4),
  mc.se    = round(c(sqrt(pce_hat * (1 - pce_hat) / nsim),
                     sqrt(fwer_hat * (1 - fwer_hat) / nsim), NA), 4),
  row.names = NULL
)

