## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(FastSurvival)
has_rpsftm <- requireNamespace("rpsftm", quietly = TRUE) &&
  requireNamespace("survival", quietly = TRUE)

## ----design-------------------------------------------------------------------
nsim <- 1000
seed <- 20261007

df <- simdata_fast(
  nsim       = nsim,
  n          = c(300, 300),
  a.time     = c(0, 18),
  a.rate     = 600 / 18,
  h01.median = list(8, 12),    # progression
  h02.median = list(30, 36),   # death without progression
  h12.median = list(12, 12),   # death after progression
  d.hazard   = -log(1 - 0.05) / 12,
  seed       = seed
)

# Single-endpoint view of PFS (k = 1) or OS (k = 2) for analysis_fast().
ep <- function(d, k) {
  data.frame(
    sim          = d$sim,
    group        = d$group,
    accrual_time = d$accrual_time,
    tte          = d[[paste0("e", k, "_tte")]],
    event        = d[[paste0("e", k, "_event")]]
  )
}

## ----pfs----------------------------------------------------------------------
pfs_cut <- cutoff_fast(df, event.looks = 380,
                       tte.col = "e1_tte", event.col = "e1_event")
pfs_res <- analysis_fast(ep(df, 1), control = 1, cutoff.looks = pfs_cut,
                         side = 1)
pfs_pos <- (pfs_res$reached & pfs_res$logrank.p <= 0.025) %in% TRUE

c(PFS_power = mean(pfs_pos),
  mean_PFS_analysis_month = mean(pfs_cut[, 1], na.rm = TRUE))

## ----switching----------------------------------------------------------------
sw_gated <- switch_fast(df, group = 1, when = "later", cutoff = pfs_cut,
                        sims = pfs_pos, aft.factor = 1.3)
sw_all   <- switch_fast(df, group = 1, when = "later", cutoff = pfs_cut,
                        aft.factor = 1.3)
sw_prog  <- switch_fast(df, group = 1, prob = 0.5, when = "intermediate",
                        median = 18)

## ----invariance---------------------------------------------------------------
os_ia0 <- analysis_fast(ep(df, 2), control = 1, cutoff.looks = pfs_cut,
                        side = 1)
os_ia1 <- analysis_fast(ep(sw_gated, 2), control = 1, cutoff.looks = pfs_cut,
                        side = 1)
all.equal(os_ia0, os_ia1)

## ----os-final-----------------------------------------------------------------
os_final <- function(d) {
  cut <- cutoff_fast(d, event.looks = 350,
                     tte.col = "e2_tte", event.col = "e2_event")
  analysis_fast(ep(d, 2), control = 1, cutoff.looks = cut,
                stat = c("logrank", "coxph"), side = 1)
}

scen <- list(
  "no switching"                   = df,
  "crossover after positive PFS"   = sw_gated,
  "crossover after PFS, all trials" = sw_all,
  "50% switch at progression"      = sw_prog
)
tab <- do.call(rbind, lapply(names(scen), function(nm) {
  d <- scen[[nm]]
  r <- os_final(d)
  data.frame(
    scenario         = nm,
    control_switched = round(mean(d$switched[d$group == 1] == 1), 3),
    OS_power         = round(mean(r$logrank.p <= 0.025, na.rm = TRUE), 3),
    mean_HR          = round(exp(mean(r$cox.coef, na.rm = TRUE)), 3),
    mean_OS_month    = round(mean(r$cutoff, na.rm = TRUE), 1)
  )
}))
tab

## ----timing-------------------------------------------------------------------
system.time({
  d0  <- simdata_fast(nsim = nsim, n = c(300, 300), a.time = c(0, 18),
                      a.rate = 600 / 18,
                      h01.median = list(8, 12), h02.median = list(30, 36),
                      h12.median = list(12, 12),
                      d.hazard = -log(1 - 0.05) / 12, seed = seed)
  pc  <- cutoff_fast(d0, event.looks = 380,
                     tte.col = "e1_tte", event.col = "e1_event")
  pr  <- analysis_fast(ep(d0, 1), control = 1, cutoff.looks = pc, side = 1)
  pos <- (pr$reached & pr$logrank.p <= 0.025) %in% TRUE
  d1  <- switch_fast(d0, group = 1, when = "later", cutoff = pc, sims = pos,
                     aft.factor = 1.3)
  r1  <- os_final(d1)
})

## ----rpsft, eval = has_rpsftm-------------------------------------------------
f   <- 1.5
h01 <- log(2) / 8
h02 <- log(2) / 30
h12 <- log(2) / 12
n_rep <- 10
big <- simdata_fast(
  nsim = n_rep, n = c(2000, 2000), a.time = c(0, 12), a.rate = 4000 / 12,
  h01.hazard = list(h01, h01 / f), h02.hazard = list(h02, h02 / f),
  h12.hazard = list(h12, h12 / f), seed = 7
)
big <- switch_fast(big, group = 1, prob = 0.6, when = "intermediate",
                   aft.factor = f, seed = 8)

suppressWarnings(suppressPackageStartupMessages({
  library(survival)
  library(rpsftm)
}))
psi_hat <- vapply(seq_len(n_rep), function(s) {
  one <- big[big$sim == s, ]
  rp_dat <- data.frame(
    time   = one$e2_surv_time,
    status = 1,
    arm    = one$group - 1,
    cens   = 1e6
  )
  # Proportion of the observed time spent on the experimental treatment.
  rp_dat$rx <- ifelse(rp_dat$arm == 1, 1,
                      ifelse(one$switched == 1,
                             (one$e2_surv_time - one$switch_time) /
                               one$e2_surv_time, 0))
  fit <- tryCatch(
    rpsftm(Surv(time, status) ~ rand(arm, rx), data = rp_dat,
           censor_time = cens),
    error = function(e) NULL
  )
  if (is.null(fit)) NA_real_ else unname(fit$psi)
}, numeric(1))

round(c(true_psi = -log(f), mean_estimate = mean(psi_hat, na.rm = TRUE),
        sd_estimate = sd(psi_hat, na.rm = TRUE),
        se_of_mean = sd(psi_hat, na.rm = TRUE) / sqrt(sum(!is.na(psi_hat)))), 4)

