## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>"
)

## ----load---------------------------------------------------------------------
library(FastSurvival)

## ----generator----------------------------------------------------------------
gen_cure_weibull <- function(nsim, n_per, accrual, shape, scale, cure,
                             drop_rate, seed) {
  set.seed(seed)
  n_tot <- 2L * n_per
  N     <- nsim * n_tot
  sim   <- rep(seq_len(nsim), each = n_tot)
  group <- rep(rep(1:2, each = n_per), times = nsim)
  accrual_time <- stats::runif(N, 0, accrual)
  surv_time    <- stats::rweibull(N, shape = shape, scale = scale[group])
  surv_time[stats::runif(N) < cure[group]] <- Inf   # cured subjects
  dropout_time <- stats::rexp(N, rate = drop_rate)
  tte <- pmin(surv_time, dropout_time)
  data.frame(sim, group, accrual_time, surv_time, dropout_time, tte,
             event = as.integer(surv_time <= dropout_time),
             calendar_time = accrual_time + tte)
}

shape <- 1.3
scale <- c(16, 18)
cure  <- c(0, 0.15)

dat <- gen_cure_weibull(nsim = 2000, n_per = 200, accrual = 12,
                        shape = shape, scale = scale, cure = cure,
                        drop_rate = 0.01, seed = 2026)
head(dat)

## ----check-generator----------------------------------------------------------
surv_model <- function(t, g) cure[g] + (1 - cure[g]) * exp(-(t / scale[g])^shape)
t_chk <- c(6, 12, 24, 36)
data.frame(
  t               = t_chk,
  control_sim     = sapply(t_chk, function(t) mean(dat$surv_time[dat$group == 1] > t)),
  control_model   = surv_model(t_chk, 1),
  treatment_sim   = sapply(t_chk, function(t) mean(dat$surv_time[dat$group == 2] > t)),
  treatment_model = surv_model(t_chk, 2)
)

## ----cutoffs------------------------------------------------------------------
cut <- cutoff_fast(dat, event.looks = c(150, 280), max.time = c(24, 48),
                   min.gap = c(NA, 9))
head(cut)
colMeans(cut)

## ----analysis-----------------------------------------------------------------
alpha_look <- c(0.002, 0.024)
set.seed(1)
res <- analysis_fast(dat, control = 1, cutoff.looks = cut,
                     stat = c("logrank", "maxcombo"), side = 1,
                     mc.alpha = alpha_look)

oc_lr <- simsummary_fast(res, p.col = "logrank.p",  alpha = alpha_look)
oc_mc <- simsummary_fast(res, p.col = "maxcombo.p", alpha = alpha_look)
data.frame(
  test             = c("Log-rank", "Max-combo"),
  reject_interim   = c(oc_lr[oc_lr$look == "1", "prob.stop.efficacy"],
                       oc_mc[oc_mc$look == "1", "prob.stop.efficacy"]),
  power            = c(oc_lr[oc_lr$look == "overall", "cum.reject"],
                       oc_mc[oc_mc$look == "overall", "cum.reject"])
)
mean(res$maxcombo.p.exact, na.rm = TRUE)

## ----switching----------------------------------------------------------------
sw <- switch_fast(dat, group = 1, when = "cutoff",
                  cutoff = cut[, 1, drop = FALSE], aft.factor = 1.5)
cut_sw <- cutoff_fast(sw, event.looks = c(150, 280), max.time = c(24, 48),
                      min.gap = c(NA, 9))
all.equal(cut_sw[, 1], cut[, 1])

res_sw <- analysis_fast(sw, control = 1, cutoff.looks = cut_sw,
                        stat = "logrank", side = 1)
all.equal(res_sw$logrank.z[res_sw$look == 1], res$logrank.z[res$look == 1])

oc_sw <- simsummary_fast(res_sw, p.col = "logrank.p", alpha = alpha_look)
c(log_rank_power_without_switching = oc_lr[oc_lr$look == "overall", "cum.reject"],
  log_rank_power_with_crossover    = oc_sw[oc_sw$look == "overall", "cum.reject"])

