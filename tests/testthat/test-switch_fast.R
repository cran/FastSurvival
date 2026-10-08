make_sw_data <- function(nsim = 20, seed = 21, d.hazard = 0.01) {
  simdata_fast(nsim = nsim, n = c(100, 100), a.time = c(0, 12),
               a.rate = 200 / 12, e.median = list(12, 18),
               d.hazard = d.hazard, seed = seed)
}

make_id_data <- function(nsim = 100, seed = 31, switch.prop = NULL,
                         h12.switch.hazard = NULL) {
  simdata_fast(nsim = nsim, n = c(200, 200), a.time = c(0, 12),
               a.rate = 400 / 12,
               h01.hazard = list(0.10, 0.07), h02.hazard = list(0.03, 0.02),
               h12.hazard = list(0.08, 0.08), switch.prop = switch.prop,
               h12.switch.hazard = h12.switch.hazard, seed = seed)
}

test_that("switch_fast: crossover at a cutoff rescales only the post-switch time", {
  df  <- make_sw_data()
  ia  <- cutoff_fast(df, event.looks = 80)
  dfs <- switch_fast(df, group = 1, when = "cutoff", cutoff = ia,
                     aft.factor = 1.5)

  open  <- unname(ia[df$sim, 1])
  s_all <- pmax(open - df$accrual_time, 0)
  # Eligible: still at risk and on study at the opening, judged on the
  # calendar scale (accrual + time), as in the analysis.
  elig  <- df$group == 1 & df$surv_time > s_all & df$dropout_time > s_all &
    df$accrual_time + df$surv_time > open &
    df$accrual_time + df$dropout_time > open
  expect_true(any(elig))
  expect_equal(dfs$switched, as.integer(elig))
  expect_equal(dfs$switch_time[elig], unname(s_all[elig]))
  expect_true(all(is.na(dfs$switch_time[!elig])))
  expect_equal(dfs$surv_time[elig],
               unname(s_all[elig] + 1.5 * (df$surv_time[elig] - s_all[elig])))
  expect_identical(dfs$surv_time[!elig], df$surv_time[!elig])
  # Observed columns are recomputed from the latent times.
  expect_equal(dfs$tte, pmin(dfs$surv_time, dfs$dropout_time))
  expect_equal(dfs$event, as.integer(dfs$surv_time <= dfs$dropout_time))
  expect_equal(dfs$calendar_time, dfs$accrual_time + dfs$tte)
  # Unchanged columns.
  expect_identical(dfs$accrual_time, df$accrual_time)
  expect_identical(dfs$dropout_time, df$dropout_time)
  expect_identical(dfs$group, df$group)
})

test_that("switch_fast: an analysis at the opening time is unchanged", {
  df  <- make_sw_data()
  ia  <- cutoff_fast(df, event.looks = 80)
  dfs <- switch_fast(df, group = 1, when = "cutoff", cutoff = ia,
                     aft.factor = 2)
  st  <- c("logrank", "coxph", "rmst")
  r0  <- analysis_fast(df,  control = 1, cutoff.looks = ia, stat = st, tau = 6)
  r1  <- analysis_fast(dfs, control = 1, cutoff.looks = ia, stat = st, tau = 6)
  expect_equal(r1, r0)
  # A later analysis is changed.
  f0 <- analysis_fast(df,  control = 1, time.looks = 60)
  f1 <- analysis_fast(dfs, control = 1, time.looks = 60)
  expect_false(isTRUE(all.equal(f0$logrank.z, f1$logrank.z)))
})

test_that("switch_fast: AFT crossover reproduces the analytic survival function", {
  nsim <- 100
  lam  <- log(2) / 12
  df   <- simdata_fast(nsim = nsim, n = c(200, 200), a.time = c(0, 6),
                       a.prop = 1, e.hazard = list(lam, lam / 1.5), seed = 41)
  # Crossover opens at calendar month 12 in every simulated trial.
  dfs  <- switch_fast(df, group = 1, when = "cutoff",
                      cutoff = rep(12, nsim), aft.factor = 2)
  ctl  <- dfs$group == 1
  s    <- 12 - dfs$accrual_time[ctl]
  # P(T > 24 | accrual) under the switch at s (independently checked).
  surv <- ifelse(24 <= s, exp(-lam * 24), exp(-lam * s - lam * (24 - s) / 2))
  # Absolute Monte Carlo tolerance (about 4 standard errors).
  expect_lt(abs(mean(dfs$surv_time[ctl] > 24) - mean(surv)), 0.015)
  # The treatment group is untouched.
  expect_identical(dfs$surv_time[!ctl], df$surv_time[!ctl])
})

test_that("switch_fast: redrawn post-switch times follow the new hazard", {
  nsim <- 100
  df   <- simdata_fast(nsim = nsim, n = c(200, 200), a.time = c(0, 6),
                       a.prop = 1, e.median = list(12, 18), seed = 42)
  sw1  <- switch_fast(df, group = 1, when = "cutoff", cutoff = rep(6, nsim),
                      median = 30)
  r1   <- (sw1$surv_time - sw1$switch_time)[sw1$switched == 1]
  expect_gt(length(r1), 5000)
  expect_lt(abs(mean(r1 > 30) - 0.5), 0.02)

  sw2 <- switch_fast(df, group = 1, when = "cutoff", cutoff = rep(6, nsim),
                     hazard = c(0.1, 0.02), time = c(0, 5, Inf))
  r2  <- (sw2$surv_time - sw2$switch_time)[sw2$switched == 1]
  expect_lt(abs(mean(r2 > 5) - exp(-0.5)), 0.02)
  expect_lt(abs(mean(r2 > 10) - exp(-0.6)), 0.02)
})

test_that("switch_fast: prob selects a random fraction of the eligible subjects", {
  nsim <- 100
  df   <- simdata_fast(nsim = nsim, n = c(200, 200), a.time = c(0, 6),
                       a.prop = 1, e.median = list(12, 18), seed = 43)
  all1 <- switch_fast(df, group = 1, when = "cutoff", cutoff = rep(6, nsim),
                      aft.factor = 1.5)
  part <- switch_fast(df, group = 1, when = "cutoff", cutoff = rep(6, nsim),
                      aft.factor = 1.5, prob = 0.3, seed = 7)
  elig <- all1$switched == 1
  expect_true(all(part$switched[!elig] == 0))
  expect_lt(abs(mean(part$switched[elig]) - 0.3), 0.02)
  # Reproducible from the seed.
  again <- switch_fast(df, group = 1, when = "cutoff", cutoff = rep(6, nsim),
                       aft.factor = 1.5, prob = 0.3, seed = 7)
  expect_identical(again, part)
})

test_that("switch_fast: switching at progression agrees with simdata_fast switch.prop", {
  idA <- make_id_data(seed = 31)
  swA <- switch_fast(idA, group = 1, prob = 0.4, when = "intermediate",
                     hazard = 0.04)
  idB <- make_id_data(seed = 32, switch.prop = list(0.4, 0),
                      h12.switch.hazard = list(0.04, 0.04))
  ctlA <- swA$group == 1 & swA$intermediate == 1
  ctlB <- idB$group == 1 & idB$intermediate == 1
  expect_lt(abs(mean(swA$switched[ctlA]) - 0.4), 0.02)
  expect_lt(abs(mean(idB$switched[ctlB]) - 0.4), 0.02)
  expect_equal(swA$switch_time[swA$switched == 1],
               swA$e1_surv_time[swA$switched == 1])
  # Same model, independent draws: the overall survival distributions agree.
  for (t in c(12, 24, 36)) {
    pA <- mean(swA$e2_surv_time[swA$group == 1] > t)
    pB <- mean(idB$e2_surv_time[idB$group == 1] > t)
    expect_lt(abs(pA - pB), 0.02)
  }
  # The first endpoint and the treatment group are untouched.
  expect_identical(swA$e1_surv_time, idA$e1_surv_time)
  expect_identical(swA$e2_surv_time[swA$group == 2],
                   idA$e2_surv_time[idA$group == 2])
})

test_that("switch_fast: crossover after a cutoff can depend on the interim result", {
  idA  <- make_id_data(nsim = 40, seed = 33)
  ia   <- cutoff_fast(idA, event.looks = 150, tte.col = "e1_tte",
                      event.col = "e1_event")
  gate <- rep(c(TRUE, FALSE), length.out = 40)
  swL  <- switch_fast(idA, group = 1, when = "later", cutoff = ia, sims = gate,
                      aft.factor = 1.4)

  off <- !gate[idA$sim]
  expect_equal(swL[off, ], idA[off, ])
  on  <- swL$switched == 1
  expect_true(any(on))
  expect_true(all(gate[swL$sim[on]]))
  open_r <- pmax(ia[swL$sim[on], 1] - swL$accrual_time[on], 0)
  expect_equal(swL$switch_time[on],
               unname(pmax(swL$e1_surv_time[on], open_r)))
  expect_identical(swL$e1_surv_time, idA$e1_surv_time)

  # Overall survival analyzed at the interim cutoff is unchanged.
  os_dat <- function(d) {
    data.frame(sim = d$sim, group = d$group, accrual_time = d$accrual_time,
               tte = d$e2_tte, event = d$e2_event)
  }
  r0 <- analysis_fast(os_dat(idA), control = 1, cutoff.looks = ia)
  r1 <- analysis_fast(os_dat(swL), control = 1, cutoff.looks = ia)
  expect_equal(r1, r0)
})

test_that("switch_fast: illness-death crossover at a cutoff accelerates both events", {
  idA <- make_id_data(nsim = 20, seed = 34)
  sw  <- switch_fast(idA, group = 1, when = "cutoff", cutoff = rep(10, 20),
                     aft.factor = 1.5)
  on  <- sw$switched == 1
  s   <- sw$switch_time[on]
  e1_old <- idA$e1_surv_time[on]
  a_on   <- idA$accrual_time[on]
  e1_exp <- ifelse(a_on + e1_old > 10, s + 1.5 * (e1_old - s), e1_old)
  expect_equal(sw$e1_surv_time[on], e1_exp)
  expect_true(all(sw$e1_surv_time <= sw$e2_surv_time))
  expect_error(switch_fast(idA, group = 1, when = "cutoff",
                           cutoff = rep(10, 20), hazard = 0.05),
               "only 'aft.factor'")
})

test_that("switch_fast: input validation", {
  df <- make_sw_data(nsim = 3)
  expect_error(switch_fast(df, group = 1, when = "intermediate",
                           aft.factor = 1.5), "needs an intermediate event")
  expect_error(switch_fast(df, group = 9, when = "cutoff", cutoff = 1:3,
                           aft.factor = 1.5), "group labels")
  expect_error(switch_fast(df, group = 1, when = "cutoff", cutoff = 1:3),
               "exactly one of")
  expect_error(switch_fast(df, group = 1, when = "cutoff", cutoff = 1:3,
                           aft.factor = 1.5, median = 10), "exactly one of")
  expect_error(switch_fast(df, group = 1, when = "cutoff", aft.factor = 1.5),
               "'cutoff' must be supplied")
  expect_error(switch_fast(df, group = 1, when = "cutoff", cutoff = 1:2,
                           aft.factor = 1.5), "one element per")
  expect_error(switch_fast(df, group = 1, when = "cutoff", cutoff = 1:3,
                           hazard = c(0.1, 0.2), time = c(0, 5)), "'time'")
  expect_error(switch_fast(df, group = 1, when = "cutoff", cutoff = 1:3,
                           hazard = 0), "positive last")
  expect_error(switch_fast(df, group = 1, when = "cutoff", cutoff = 1:3,
                           hazard = 0.1, time = c(0, Inf)), "piecewise")
  expect_error(switch_fast(df, group = 1, when = "cutoff", cutoff = 1:3,
                           aft.factor = 1.5, sims = TRUE), "'sims'")
  expect_error(switch_fast(df, group = 1, when = "cutoff", cutoff = 1:3,
                           aft.factor = 1.5, stream = 1), "requires 'seed'")
  id <- make_id_data(nsim = 3)
  expect_error(switch_fast(id, group = 1, when = "intermediate",
                           aft.factor = 1.5, sims = c(TRUE, FALSE, TRUE)),
               "'sims' requires")
  expect_error(switch_fast(id, group = 1, when = "intermediate",
                           cutoff = 1:3, aft.factor = 1.5),
               "'cutoff' is used only")
  expect_error(switch_fast(data.frame(sim = 1, group = 1, accrual_time = 0),
                           group = 1, when = "cutoff", cutoff = 1,
                           aft.factor = 2), "latent columns")
})

test_that("switch_fast: a named cutoff vector is matched by name", {
  df  <- make_sw_data(nsim = 5)
  ia  <- cutoff_fast(df, event.looks = 60)
  ref <- switch_fast(df, group = 1, when = "cutoff", cutoff = ia,
                     aft.factor = 1.5)
  v <- ia[, 1]
  expect_identical(names(v), rownames(ia))
  expect_identical(switch_fast(df, group = 1, when = "cutoff", cutoff = rev(v),
                               aft.factor = 1.5), ref)
  expect_error(switch_fast(df, group = 1, when = "cutoff",
                           cutoff = stats::setNames(v, paste0("x", 1:5)),
                           aft.factor = 1.5), "names")
})

test_that("switch_fast: subjects without a finite event or dropout time stay censored", {
  df <- simdata_fast(nsim = 5, n = c(100, 100), a.time = c(0, 12),
                     a.rate = 200 / 12, e.hazard = list(c(0.1, 0), c(0.1, 0)),
                     e.time = c(0, 24, Inf), seed = 6)
  sw <- switch_fast(df, group = 1, when = "cutoff", cutoff = rep(18, 5),
                    aft.factor = 1.5)
  expect_true(any(sw$switched == 1 & !is.finite(sw$surv_time)))
  expect_equal(sw$event, as.integer(is.finite(sw$surv_time)))
})
