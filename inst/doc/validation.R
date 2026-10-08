## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>"
)
have_survival <- requireNamespace("survival", quietly = TRUE)
have_survRM2  <- requireNamespace("survRM2", quietly = TRUE)
have_survAH   <- requireNamespace("survAH", quietly = TRUE)
have_nph      <- requireNamespace("nph", quietly = TRUE)
have_nphRCT   <- requireNamespace("nphRCT", quietly = TRUE)
have_simtrial <- requireNamespace("simtrial", quietly = TRUE)

## ----load---------------------------------------------------------------------
library(FastSurvival)

## ----data, eval = have_survival-----------------------------------------------
library(survival)

# German Breast Cancer Study Group cohort
str(gbsg[, c("rfstime", "status", "hormon")])
table(hormon = gbsg$hormon)

## ----km, eval = have_survival-------------------------------------------------
# survfit_fast assumes time-sorted input; the C++ core groups tied times in the
# survival::survfit risk-set convention, so the data only need ordering by time.
ord  <- order(gbsg$rfstime)
fast <- survfit_fast(gbsg$rfstime[ord], gbsg$status[ord],
                     t_eval = 1000, conf.type = "log-log")

fit  <- survfit(Surv(rfstime, status) ~ 1, data = gbsg)
ref  <- summary(fit, times = 1000)

data.frame(
  quantity = c("survival", "std.err"),
  fast     = c(unclass(fast)["surv"], unclass(fast)["std.err"]),
  survival = c(ref$surv, ref$std.err),
  row.names = NULL
)

## ----logrank, eval = have_survival--------------------------------------------
fast_lr <- survdiff_fast(gbsg$rfstime, gbsg$status, gbsg$hormon,
                         control = 0, side = 1)

ref_lr  <- survdiff(Surv(rfstime, status) ~ hormon, data = gbsg)

c(fast = as.numeric(fast_lr)^2, survival = ref_lr$chisq)

## ----wlr, eval = have_survival && have_nph------------------------------------
fh_grid <- data.frame(rho = c(0, 1, 0, 1), gamma = c(1, 0, 0, 1))

do.call(rbind, lapply(seq_len(nrow(fh_grid)), function(i) {
  r <- fh_grid$rho[i]
  g <- fh_grid$gamma[i]
  fast <- as.numeric(survdiff_fast(gbsg$rfstime, gbsg$status, gbsg$hormon,
                                   control = 0, side = 1,
                                   weight = "fh", rho = r, gamma = g))
  nph_chisq <- nph::logrank.test(gbsg$rfstime, gbsg$status, gbsg$hormon,
                                 rho = r, gamma = g)$test$Chisq
  data.frame(rho = r, gamma = g, fast = fast^2, nph = nph_chisq)
}))

## ----cox, eval = have_survival------------------------------------------------
fast_hr <- coxph_fast(gbsg$rfstime, gbsg$status, gbsg$hormon,
                      control = 0, side = 1)

ref_cox <- coxph(Surv(rfstime, status) ~ hormon, data = gbsg,
                 ties = "breslow")

c(fast = unclass(fast_hr)["coef"], cox = unname(coef(ref_cox)))
cox_diff <- abs(unname(unclass(fast_hr)["coef"]) - unname(coef(ref_cox)))

## ----cox-strata, eval = have_survival-----------------------------------------
fast_hr_s <- coxph_fast(gbsg$rfstime, gbsg$status, gbsg$hormon,
                        control = 0, strata = gbsg$grade)

ref_cox_s <- coxph(Surv(rfstime, status) ~ hormon + strata(grade),
                   data = gbsg, ties = "breslow")

c(fast = unclass(fast_hr_s)["coef"], cox = unname(coef(ref_cox_s)))

## ----rmst, eval = have_survival && have_survRM2-------------------------------
library(survRM2)

tau <- 1000

fast_rmst <- rmst_fast(gbsg$rfstime, gbsg$status, gbsg$hormon,
                       control = 0, tau = tau, side = 1)

ref_rmst  <- rmst2(time = gbsg$rfstime, status = gbsg$status,
                   arm = gbsg$hormon, tau = tau)

c(fast = unclass(fast_rmst)["diff"],
  survRM2 = ref_rmst$unadjusted.result[1, 1])

## ----wmst, eval = have_survival-----------------------------------------------
tau1 <- 200
tau2 <- 1000

fast_wm <- wmst_fast(gbsg$rfstime, gbsg$status, gbsg$hormon,
                     control = 0, side = 1, tau1 = tau1, tau2 = tau2)

# Integral of the Kaplan-Meier step function over [lo, hi] for one group.
km_window_area <- function(gi, lo, hi) {
  sf <- survfit(Surv(rfstime, status) ~ 1, data = gbsg[gbsg$hormon == gi, ])
  tk <- c(0, sf$time)
  sk <- c(1, sf$surv)
  brk  <- sort(unique(c(lo, hi, sf$time[sf$time > lo & sf$time < hi])))
  area <- 0
  for (b in seq_len(length(brk) - 1L)) {
    u  <- brk[b]
    su <- sk[max(which(tk <= u))]
    area <- area + su * (brk[b + 1L] - u)
  }
  area
}

data.frame(
  group    = c("control", "treatment"),
  fast     = unclass(fast_wm)[c("wmst.control", "wmst.treatment")],
  survfit  = c(km_window_area(0, tau1, tau2), km_window_area(1, tau1, tau2)),
  row.names = NULL
)

## ----wkm, eval = have_survival------------------------------------------------
tmax <- max(gbsg$rfstime)

fast_wk <- wkm_fast(gbsg$rfstime, gbsg$status, gbsg$hormon,
                    control = 0, side = 1, weight = "constant")

# tau = tmax extends past the shorter group's follow-up on purpose, to cover
# the whole observed range as the constant-weight test does; rmst_fast()
# warns about this, so the warning is suppressed here.
ref_rm  <- suppressWarnings(
  rmst_fast(gbsg$rfstime, gbsg$status, gbsg$hormon,
            control = 0, tau = tmax, side = 1))

c(wkm.constant = unclass(fast_wk)["wdiff"], rmst = unclass(ref_rm)["diff"])

## ----milestone, eval = have_survival------------------------------------------
tstar <- 1000

fast_ms <- milestone_fast(gbsg$rfstime, gbsg$status, gbsg$hormon,
                          control = 0, tau = tstar, method = "loglog",
                          side = 1)

fit_g <- survfit(Surv(rfstime, status) ~ hormon, data = gbsg)
ref_g <- summary(fit_g, times = tstar)

data.frame(
  group    = c("control", "treatment"),
  fast     = fast_ms$surv[c("control", "treatment")],
  survival = ref_g$surv,
  row.names = NULL
)

## ----medsurv, eval = have_survival--------------------------------------------
fast_med <- medsurv_fast(gbsg$rfstime, gbsg$status, gbsg$hormon,
                         control = 0, side = 1, method = "nph")

fit_med <- survfit(Surv(rfstime, status) ~ hormon, data = gbsg)
med_ref <- summary(fit_med)$table[, "median"]

data.frame(
  group    = c("control", "treatment"),
  fast     = unclass(fast_med)[c("median.control", "median.treatment")],
  survival = c(med_ref[1], med_ref[2]),
  row.names = NULL
)

## ----ahsw, eval = have_survival && have_survAH--------------------------------
library(survAH)

tau <- 1000

fast_ah <- ahsw_fast(gbsg$rfstime, gbsg$status, gbsg$hormon,
                     control = 0, tau = tau, side = 1)

ref_ah  <- ah2(time = gbsg$rfstime, status = gbsg$status,
               arm = gbsg$hormon, tau = tau)

data.frame(
  quantity = c("AH (control)", "AH (treatment)"),
  fast     = unclass(fast_ah)[c("ah.ctrl", "ah.trt")],
  survAH   = c(ref_ah$ah["AH (arm0)", "Est."],
               ref_ah$ah["AH (arm1)", "Est."]),
  row.names = NULL
)

## ----ahr, eval = have_survival------------------------------------------------
tau <- 1000

fast_ahr <- ahr_fast(gbsg$rfstime, gbsg$status, gbsg$hormon,
                     control = 0, tau = tau, side = 1)

g  <- gbsg$hormon
ev <- c(gbsg$rfstime[g == 0 & gbsg$status == 1],
        gbsg$rfstime[g == 1 & gbsg$status == 1])
grid <- sort(unique(c(0, ev[ev <= tau], tau)))

km_on_grid <- function(gi) {
  sf <- survfit(Surv(rfstime, status) ~ 1, data = gbsg[g == gi, ])
  approxfun(sf$time, sf$surv, method = "constant",
            yleft = 1, rule = 2, f = 0)(grid)
}

S0  <- km_on_grid(0)
S1  <- km_on_grid(1)
m   <- length(grid)
dS0 <- S0 - c(1, S0[-m])
GL  <- S0[m] * S1[m]
ref_theta_ctrl <- -sum(S1 * dS0) / (1 - GL)

data.frame(
  quantity  = c("theta (control)", "theta (treatment)", "AHR"),
  fast      = c(fast_ahr$theta[[1]], fast_ahr$theta[[2]], fast_ahr$ahr),
  reference = c(ref_theta_ctrl, 1 - ref_theta_ctrl,
                (1 - ref_theta_ctrl) / ref_theta_ctrl),
  row.names = NULL
)

## ----maxcombo, eval = have_survival && have_nph-------------------------------
rho   <- c(0, 0, 1)
gamma <- c(0, 1, 0)

fast_mc <- maxcombo_fast(gbsg$rfstime, gbsg$status, gbsg$hormon,
                         control = 0, side = 1, rho = rho, gamma = gamma)

nph_mc <- nph::logrank.maxtest(gbsg$rfstime, gbsg$status, gbsg$hormon)

data.frame(
  weight = c("FH(0,0)", "FH(0,1)", "FH(1,0)"),
  fast   = abs(attr(fast_mc, "z")),
  nph    = abs(nph_mc$tests$z),
  row.names = NULL
)

## ----rmw, eval = have_survival && have_nphRCT---------------------------------
gbsg_df <- data.frame(
  time  = gbsg$rfstime,
  event = gbsg$status,
  arm   = factor(ifelse(gbsg$hormon == 0, "control", "experimental"),
                 levels = c("control", "experimental"))
)

fit_rmw <- rmw_fast(gbsg_df$time, gbsg_df$event, gbsg_df$arm,
                    control = "control", side = 1, s_star = 0.5)

z_lr_nph <- nphRCT::wlrt(Surv(time, event) ~ arm, data = gbsg_df,
                         method = "mw", s_star = 1)$z
z_mw_nph <- nphRCT::wlrt(Surv(time, event) ~ arm, data = gbsg_df,
                         method = "mw", s_star = 0.5)$z

data.frame(
  component = c("log-rank (s_star = 1)", "modestly-weighted (s_star = 0.5)"),
  fast      = abs(attr(fit_rmw, "z")),
  nphRCT    = abs(c(z_lr_nph, z_mw_nph)),
  row.names = NULL
)

## ----cutoff, eval = have_simtrial---------------------------------------------
sim_c <- simdata_fast(nsim = 50, n = c(150, 150), a.time = c(0, 12),
                      a.rate = 300 / 12, e.median = list(12, 16),
                      d.hazard = 0.01, seed = 2026)
cut_c <- cutoff_fast(sim_c, event.looks = 150, time.looks = 18,
                     max.time = 30, min.enrolled = 250, min.followup = 6)
ref_c <- vapply(seq_len(50), function(s) {
  d <- sim_c[sim_c$sim == s, ]
  x <- data.frame(enroll_time = d$accrual_time,
                  cte = d$accrual_time + d$tte,
                  fail = d$event, stratum = "All")
  simtrial::get_analysis_date(x, planned_calendar_time = 18,
                              target_event_overall = 150,
                              max_extension_for_target_event = 30,
                              min_n_overall = 250, min_followup = 6)
}, numeric(1))
c(max_abs_difference = max(abs(cut_c[, 1] - ref_c)))

## ----switching----------------------------------------------------------------
lam   <- log(2) / 12
sim_s <- simdata_fast(nsim = 200, n = c(250, 250), a.time = c(0, 6),
                      a.prop = 1, e.hazard = list(lam, lam / 1.5),
                      seed = 2027)
sw    <- switch_fast(sim_s, group = 1, when = "cutoff",
                     cutoff = rep(12, 200), aft.factor = 2)
ctl   <- sw$group == 1
s_sw  <- 12 - sw$accrual_time[ctl]
t_grid <- c(6, 12, 18, 24, 36)
data.frame(
  t         = t_grid,
  empirical = round(sapply(t_grid, function(t) mean(sw$surv_time[ctl] > t)), 4),
  analytic  = round(sapply(t_grid, function(t) {
    mean(ifelse(t <= s_sw, exp(-lam * t),
                exp(-lam * s_sw - lam * (t - s_sw) / 2)))
  }), 4)
)

