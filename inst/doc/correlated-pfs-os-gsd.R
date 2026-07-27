## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6.5,
  fig.height = 4.2,
  fig.align = "center"
)
library(FastSurvival)
has_gsd <- requireNamespace("gsDesign", quietly = TRUE)

## ----design-------------------------------------------------------------------
nsim <- 5000
seed <- 20260611

# Per-group sample sizes (group 1 = control, group 2 = treatment)
n_control <- 300
n_treat   <- 300
r_alloc   <- n_treat / n_control          # allocation ratio treatment:control

# Piecewise-uniform accrual: 20 per month, then 40 per month (600 enrolled by month 18)
a_time <- c(0, 6)
a_rate <- c(20, 40)

# Median PFS and OS by group, in months
mst_PFS_C <- 6
mst_PFS_T <- 9
mst_OS_C  <- 14
mst_OS_T  <- 18

# Exponential hazards implied by the medians
lam_PFS_C <- log(2) / mst_PFS_C
lam_PFS_T <- log(2) / mst_PFS_T
lam_OS_C  <- log(2) / mst_OS_C
lam_OS_T  <- log(2) / mst_OS_T

# Fleischer maximal-independence: PFS hazard = TTP hazard + OS hazard,
# so the latent TTP hazard is the PFS hazard minus the OS hazard.
lam_TTP_C <- lam_PFS_C - lam_OS_C
lam_TTP_T <- lam_PFS_T - lam_OS_T
stopifnot(lam_TTP_C > 0, lam_TTP_T > 0)

# Non-informative exponential dropout, 10 percent per year on the month scale
eta <- -log(1 - 0.10) / 12

# Event-driven looks (cumulative PFS event counts):
#   look 1 = futility, look 2 = efficacy interim, look 3 = efficacy final
d_PFS <- c(200, 300, 400)
L <- length(d_PFS)

# One-sided overall significance level
alpha_total <- 0.025

c(HR_PFS = lam_PFS_T / lam_PFS_C, HR_OS = lam_OS_T / lam_OS_C)

## ----generate-----------------------------------------------------------------
# A single call returns both correlated endpoints. e1 is PFS (the state-0 exit:
# progression or death), e2 is OS (the terminal event), and dropout_time is the
# shared non-informative dropout. The observed columns e1_tte / e1_event and
# e2_tte / e2_event are each already censored at the shared dropout time.
df <- simdata_fast(
  nsim       = nsim,
  n          = c(n_control, n_treat),
  a.time     = a_time,
  a.rate     = a_rate,
  h01.hazard = list(lam_TTP_C, lam_TTP_T),
  h02.hazard = list(lam_OS_C, lam_OS_T),
  d.hazard   = eta,
  seed       = seed
)

# Observed PFS endpoint (e1).
pfs_dat <- data.frame(
  sim          = df$sim,
  group        = df$group,
  accrual_time = df$accrual_time,
  tte          = df$e1_tte,
  event        = df$e1_event
)

# Observed OS endpoint (e2).
os_dat <- data.frame(
  sim          = df$sim,
  group        = df$group,
  accrual_time = df$accrual_time,
  os_tte       = df$e2_tte,
  os_event     = df$e2_event
)

## ----pfs-analysis-------------------------------------------------------------
pfs_res <- analysis_fast(
  pfs_dat,
  control     = 1,
  event.looks = d_PFS,
  stat        = "logrank",
  side        = 1
)

# Per-simulation calendar cutoffs A_l (rows = simulation, columns = look)
A_mat <- matrix(pfs_res$cutoff, nrow = nsim, ncol = L, byrow = TRUE)

# With this sample size and these event targets the looks are reached in
# essentially every simulation; warn if any target is not reached.
if (!all(pfs_res$reached)) {
  warning("Some PFS event targets were not reached; increase n or lower d_PFS.")
}

# Standardized PFS log-rank statistics (natural sign: benefit is negative)
Z_PFS <- matrix(pfs_res$logrank.z, nrow = nsim, ncol = L, byrow = TRUE)

## ----os-analysis--------------------------------------------------------------
big <- max(A_mat, na.rm = TRUE) + 1

os_blocks <- vector("list", L)
for (l in seq_len(L)) {
  A_l    <- A_mat[os_dat$sim, l]                        # per-row cutoff at look l
  follow <- pmax(0, A_l - os_dat$accrual_time)          # administrative follow-up
  obs_t  <- pmin(os_dat$os_tte, follow)                 # observed OS time at the look
  obs_e  <- os_dat$os_event *
    as.integer(os_dat$accrual_time + os_dat$os_tte <= A_l)

  cut_l <- data.frame(
    sim          = os_dat$sim,
    group        = os_dat$group,
    accrual_time = 0,
    tte          = obs_t,
    event        = obs_e
  )
  res_l <- analysis_fast(
    cut_l, control = 1, time.looks = big, stat = "logrank", side = 1
  )

  # Relabel as look l and report the PFS-driven cutoff for the timing summary;
  # keep only the columns that are meaningful for the censored OS data.
  res_l$look       <- l
  res_l$look.value <- d_PFS[l]
  res_l$cutoff     <- A_mat[, l]
  os_blocks[[l]] <- res_l[, c("sim", "look", "look.value", "cutoff",
                              "n.event", "logrank.z", "logrank.chisq",
                              "logrank.p")]
}

os_res <- do.call(rbind, os_blocks)
os_res <- os_res[order(os_res$sim, os_res$look), ]

# Standardized OS log-rank statistics
Z_OS <- matrix(NA_real_, nsim, L)
for (l in seq_len(L)) Z_OS[, l] <- os_blocks[[l]]$logrank.z

## ----event-counts-------------------------------------------------------------
d_PFS_emp <- colMeans(matrix(pfs_res$n.event, nrow = nsim, ncol = L, byrow = TRUE))
d_OS_emp  <- colMeans(matrix(os_res$n.event,  nrow = nsim, ncol = L, byrow = TRUE))
A_mean    <- colMeans(A_mat)

data.frame(
  look          = seq_len(L),
  d_PFS_target  = d_PFS,
  d_PFS_mean    = round(d_PFS_emp),
  d_OS_mean     = round(d_OS_emp),
  cutoff_months = round(A_mean, 1)
)

## ----boundaries, eval = has_gsd-----------------------------------------------
timing_PFS_eff <- c(d_PFS[2] / d_PFS[3], 1)
timing_OS_eff  <- c(d_OS_emp[2] / d_OS_emp[3], 1)

gs_PFS <- gsDesign::gsDesign(
  k = 2, test.type = 1, alpha = alpha_total,
  timing = timing_PFS_eff, sfu = gsDesign::sfLDOF
)
gs_OS <- gsDesign::gsDesign(
  k = 2, test.type = 1, alpha = alpha_total,
  timing = timing_OS_eff, sfu = gsDesign::sfLDPocock
)

b_PFS <- gs_PFS$upper$bound      # positive z-boundaries (interim, final)
b_OS  <- gs_OS$upper$bound

# Efficacy boundaries on the natural-sign scale; NA at the futility-only look 1.
eff_PFS <- c(NA, -b_PFS[1], -b_PFS[2])
eff_OS  <- c(NA, -b_OS[1],  -b_OS[2])

# Non-binding PFS futility at look 1: stop if the observed PFS hazard ratio is at
# or above futility_HR. On the natural-sign scale this is logrank.z at or above
# log(futility_HR) * sqrt(r * d) / (1 + r).
futility_HR <- 1.0
fut1 <- log(futility_HR) * sqrt(r_alloc * d_PFS[1]) / (1 + r_alloc)
fut_PFS <- c(fut1, NA, NA)

data.frame(
  look        = seq_len(L),
  PFS_efficacy = round(eff_PFS, 3),
  PFS_futility = round(fut_PFS, 3),
  OS_efficacy  = round(eff_OS, 3)
)

## ----oc-pfs, eval = has_gsd---------------------------------------------------
oc_PFS <- simsummary_fast(
  pfs_res,
  eff.col  = "logrank.z", efficacy = eff_PFS,
  fut.col  = "logrank.z", futility = fut_PFS,
  direction = "lower"
)
oc_PFS

## ----oc-os, eval = has_gsd----------------------------------------------------
oc_OS <- simsummary_fast(
  os_res,
  eff.col = "logrank.z", efficacy = eff_OS,
  direction = "lower"
)
oc_OS

## ----joint, eval = has_gsd----------------------------------------------------
# First efficacy-crossing look for a matrix of statistics, honouring an optional
# futility rule that stops the trial without a rejection.
reject_look <- function(Z, efficacy, futility = NULL) {
  nl   <- ncol(Z)
  ns   <- nrow(Z)
  stop_flag <- logical(ns)
  out  <- rep(NA_integer_, ns)
  for (l in seq_len(nl)) {
    active <- !stop_flag
    if (!is.na(efficacy[l])) {
      hit <- active & !is.na(Z[, l]) & Z[, l] <= efficacy[l]
      out[hit] <- l
      stop_flag[hit] <- TRUE
    }
    if (!is.null(futility) && !is.na(futility[l])) {
      fut_hit <- !stop_flag & active & !is.na(Z[, l]) & Z[, l] >= futility[l]
      stop_flag[fut_hit] <- TRUE
    }
  }
  out
}

pfs_look <- reject_look(Z_PFS, eff_PFS, fut_PFS)
os_look  <- reject_look(Z_OS,  eff_OS)

pfs_reject <- !is.na(pfs_look)
os_reject  <- !is.na(os_look)

# Hierarchical gatekeeping: OS is claimed only if PFS is claimed at the same or
# an earlier look.
joint_reject <- pfs_reject & os_reject & (os_look >= pfs_look)

power_table <- data.frame(
  quantity = c("PFS power (marginal)",
               "OS power (marginal)",
               "OS power (after PFS, hierarchical)",
               "Joint power (both endpoints)"),
  value = c(mean(pfs_reject),
            mean(os_reject),
            mean(joint_reject),
            mean(joint_reject))
)
power_table$value <- round(power_table$value, 3)
power_table

## ----correlation--------------------------------------------------------------
Zc <- cbind(Z_PFS, Z_OS)
colnames(Zc) <- c(paste0("PFS_", seq_len(L)), paste0("OS_", seq_len(L)))
cor_mat <- cor(Zc, use = "pairwise.complete.obs")
round(cor_mat, 3)

## ----cross-endpoint-----------------------------------------------------------
cross_block <- cor_mat[paste0("PFS_", seq_len(L)), paste0("OS_", seq_len(L)),
                       drop = FALSE]
data.frame(
  look = seq_len(L),
  cor_PFS_OS = round(diag(cross_block), 3)
)

## ----within-endpoint----------------------------------------------------------
pairs_idx <- rbind(c(1, 2), c(1, 3), c(2, 3))

within_tab <- function(Zmat, d_counts, tag) {
  emp <- apply(pairs_idx, 1, function(p) {
    stats::cor(Zmat[, p[1]], Zmat[, p[2]], use = "complete.obs")
  })
  theo <- apply(pairs_idx, 1, function(p) {
    sqrt(min(d_counts[p]) / max(d_counts[p]))
  })
  data.frame(
    endpoint   = tag,
    look_pair  = paste0(pairs_idx[, 1], "-", pairs_idx[, 2]),
    empirical  = round(emp, 3),
    closed_form = round(theo, 3)
  )
}

rbind(
  within_tab(Z_PFS, d_PFS,    "PFS"),
  within_tab(Z_OS,  d_OS_emp, "OS")
)

## ----scatter, fig.cap = "Standardized log-rank statistics at the final look."----
plot(
  Z_PFS[, L], Z_OS[, L],
  pch = 16, col = grDevices::adjustcolor("steelblue", alpha.f = 0.25),
  xlab = expression(Z[PFS]), ylab = expression(Z[OS]),
  main = "Final-look log-rank statistics"
)
abline(h = 0, v = 0, col = "grey70", lty = 3)

