# Tests for the illness-death (two correlated endpoints, optional switching)
# path of simdata_fast, dispatched to simdata_core_id. The single-endpoint path
# is covered by test-simdata_fast.R and is left unchanged by this feature.

id_cols <- c("sim", "group", "accrual_time", "e1_surv_time", "e2_surv_time",
             "dropout_time", "e1_tte", "e1_event", "e2_tte", "e2_event",
             "e1_calendar_time", "e2_calendar_time", "intermediate",
             "switched", "switch_time")

test_that("simdata_fast (illness-death): structure and schema, two groups", {
  dat <- simdata_fast(
    nsim = 50, n = c(100, 120), a.time = c(0, 12), a.prop = 1,
    h01.median = list(8, 12), h02.median = list(24, 30),
    d.median = list(36, 36), seed = 1
  )

  expect_s3_class(dat, "data.frame")
  expect_equal(nrow(dat), 50 * (100 + 120))
  expect_identical(names(dat), id_cols)
  expect_setequal(unique(dat$group), c(1, 2))

  # Endpoint definitions and the shared dropout time.
  expect_true(all(dat$e1_tte == pmin(dat$e1_surv_time, dat$dropout_time)))
  expect_true(all(dat$e2_tte == pmin(dat$e2_surv_time, dat$dropout_time)))
  expect_true(all(dat$e1_event == as.integer(dat$e1_surv_time <= dat$dropout_time)))
  expect_true(all(dat$e2_event == as.integer(dat$e2_surv_time <= dat$dropout_time)))
  expect_equal(dat$e1_calendar_time, dat$accrual_time + dat$e1_tte, tolerance = 1e-12)
  expect_equal(dat$e2_calendar_time, dat$accrual_time + dat$e2_tte, tolerance = 1e-12)

  # The terminal endpoint never precedes the first endpoint.
  expect_true(all(dat$e2_surv_time >= dat$e1_surv_time))

  # Indicator columns are 0/1; with no switching every switch field is empty.
  expect_true(all(dat$intermediate %in% c(0L, 1L)))
  expect_true(all(dat$e1_event %in% c(0L, 1L)))
  expect_true(all(dat$e2_event %in% c(0L, 1L)))
  expect_true(all(dat$switched == 0L))
  expect_true(all(is.na(dat$switch_time)))
})

test_that("simdata_fast (illness-death): rows are interleaved in (sim, group) order", {
  nc <- 30; nt <- 40
  dat <- simdata_fast(
    nsim = 10, n = c(nc, nt), a.time = c(0, 10), a.prop = 1,
    h01.median = list(8, 10), h02.median = list(20, 24), seed = 2
  )
  for (s in 1:10) {
    block <- dat[dat$sim == s, ]
    expect_equal(block$group, c(rep(1, nc), rep(2, nt)))
  }
})

test_that("simdata_fast (illness-death): reproducible from seed, sensitive to seed", {
  args <- list(nsim = 20, n = c(50, 50), a.time = c(0, 8), a.prop = 1,
               h01.median = list(8, 10), h02.median = list(24, 28),
               d.median = list(30, 30))
  a  <- do.call(simdata_fast, c(args, list(seed = 123)))
  b  <- do.call(simdata_fast, c(args, list(seed = 123)))
  c2 <- do.call(simdata_fast, c(args, list(seed = 124)))

  expect_equal(a$e1_surv_time, b$e1_surv_time, tolerance = 0)
  expect_equal(a$e2_surv_time, b$e2_surv_time, tolerance = 0)
  expect_equal(a$accrual_time, b$accrual_time, tolerance = 0)
  expect_false(isTRUE(all.equal(a$e2_surv_time, c2$e2_surv_time)))
})

test_that("simdata_fast (illness-death): single-group mode", {
  dat <- simdata_fast(
    nsim = 30, n = 80, a.time = c(0, 12), a.prop = 1,
    h01.hazard = 0.10, h02.hazard = 0.05, seed = 3
  )
  expect_equal(nrow(dat), 30 * 80)
  expect_setequal(unique(dat$group), 1)
  expect_true(all(is.infinite(dat$dropout_time)))   # no dropout supplied
  expect_true(all(dat$e1_event == 1L))
  expect_true(all(dat$e2_event == 1L))
})

test_that("simdata_fast (illness-death): reduces to Fleischer Theorem 1", {
  # Constant hazards, no switching, h12 defaults to h02. Then the first endpoint
  # is Exp(lam1 + lam2), the terminal endpoint is Exp(lam2), the intermediate
  # fraction is lam1 / (lam1 + lam2), and Corr(e1, e2) = median(e1) / median(e2)
  # = lam2 / (lam1 + lam2) (Fleischer 2009, Theorem 1).
  lam1 <- 0.10; lam2 <- 0.05
  dat <- simdata_fast(
    nsim = 1, n = 100000, a.time = c(0, 1), a.prop = 1,
    h01.hazard = lam1, h02.hazard = lam2, seed = 4
  )

  expect_equal(1 / mean(dat$e1_surv_time), lam1 + lam2, tolerance = 0.02)
  expect_equal(1 / mean(dat$e2_surv_time), lam2, tolerance = 0.02)
  expect_equal(mean(dat$intermediate), lam1 / (lam1 + lam2), tolerance = 0.02)
  # The sample correlation of two right-skewed (exponential) endpoints has a
  # larger sampling error than the marginal means, so it uses a looser MC
  # tolerance; the theoretical value lam2 / (lam1 + lam2) is exact.
  expect_equal(cor(dat$e1_surv_time, dat$e2_surv_time),
               lam2 / (lam1 + lam2), tolerance = 0.05)
})

test_that("simdata_fast (illness-death): post-event hazard h12 is recovered (clock-reset)", {
  # Among subjects with an intermediate event and no switching, the post-event
  # survival e2 - e1 is Exp(h12) measured from the intermediate event.
  lam12 <- 0.02
  dat <- simdata_fast(
    nsim = 1, n = 100000, a.time = c(0, 1), a.prop = 1,
    h01.hazard = 0.10, h02.hazard = 0.05, h12.hazard = lam12, seed = 5
  )
  prog <- dat$intermediate == 1L
  post <- dat$e2_surv_time[prog] - dat$e1_surv_time[prog]
  expect_gt(sum(prog), 1000)
  expect_true(all(post > 0))
  expect_equal(1 / mean(post), lam12, tolerance = 0.03)
})

test_that("simdata_fast (illness-death): treatment switching works", {
  # Control switches with probability 0.4 at the intermediate event; treatment
  # does not switch. Switchers follow a separate, more favourable post-event
  # hazard h12.switch.
  lam12_sw <- 0.015
  dat <- simdata_fast(
    nsim = 1, n = c(60000, 60000), a.time = c(0, 1), a.prop = 1,
    h01.hazard = list(0.10, 0.07), h02.hazard = list(0.05, 0.04),
    switch.prop = list(0.4, 0),
    h12.switch.hazard = list(lam12_sw, lam12_sw),
    seed = 6
  )

  # Only progressors can switch; switching implies an intermediate event.
  expect_true(all(dat$switched[dat$intermediate == 0L] == 0L))
  # No switching in the treatment group.
  expect_true(all(dat$switched[dat$group == 2L] == 0L))

  # Among control progressors, the switched fraction is about 0.4.
  ctrl_prog <- dat$group == 1L & dat$intermediate == 1L
  expect_equal(mean(dat$switched[ctrl_prog]), 0.4, tolerance = 0.02)

  # The switch time equals the intermediate-event time for switchers, NA else.
  sw <- dat$switched == 1L
  expect_equal(dat$switch_time[sw], dat$e1_surv_time[sw], tolerance = 0)
  expect_true(all(is.na(dat$switch_time[!sw])))

  # Switchers' post-event survival recovers the switch hazard (clock-reset).
  post_sw <- dat$e2_surv_time[sw] - dat$e1_surv_time[sw]
  expect_gt(sum(sw), 5000)
  expect_equal(1 / mean(post_sw), lam12_sw, tolerance = 0.03)
})

test_that("simdata_fast (illness-death): two-group marginal hazards are recovered", {
  dat <- simdata_fast(
    nsim = 1, n = c(60000, 60000), a.time = c(0, 1), a.prop = 1,
    h01.hazard = list(0.10, 0.07), h02.hazard = list(0.05, 0.04),
    seed = 7
  )
  # First endpoint is Exp(h01 + h02); terminal endpoint is Exp(h02) per group.
  e1c <- 1 / mean(dat$e1_surv_time[dat$group == 1])
  e1t <- 1 / mean(dat$e1_surv_time[dat$group == 2])
  e2c <- 1 / mean(dat$e2_surv_time[dat$group == 1])
  e2t <- 1 / mean(dat$e2_surv_time[dat$group == 2])
  expect_equal(e1c, 0.15, tolerance = 0.02)
  expect_equal(e1t, 0.11, tolerance = 0.02)
  expect_equal(e2c, 0.05, tolerance = 0.02)
  expect_equal(e2t, 0.04, tolerance = 0.02)
})

test_that("simdata_fast (illness-death): dropout censors both endpoints together", {
  dat <- simdata_fast(
    nsim = 1, n = 40000, a.time = c(0, 1), a.prop = 1,
    h01.hazard = 0.10, h02.hazard = 0.05, h12.hazard = 0.03,
    d.median = 12, seed = 8
  )
  # Some terminal events are censored by dropout.
  expect_lt(mean(dat$e2_event), 1)
  expect_gt(mean(dat$e2_event), 0.3)
  # A censored terminal endpoint is capped at the shared dropout time.
  cens <- dat$e2_surv_time > dat$dropout_time
  expect_true(all(dat$e2_event[cens] == 0L))
  expect_equal(dat$e2_tte[cens], dat$dropout_time[cens], tolerance = 0)
  # A first-endpoint event always implies the terminal endpoint is observed no
  # earlier, so e1 censoring implies e2 censoring at the same dropout time.
  expect_true(all(dat$e1_tte <= dat$e2_tte))
})

test_that("simdata_fast (illness-death): input validation", {
  # h01 supplied without h02
  expect_error(
    simdata_fast(nsim = 5, n = 50, a.time = c(0, 1), a.prop = 1,
                 h01.hazard = 0.1),
    "h02")
  # switch.prop positive without a switch hazard
  expect_error(
    simdata_fast(nsim = 5, n = c(50, 50), a.time = c(0, 1), a.prop = 1,
                 h01.hazard = list(0.1, 0.08), h02.hazard = list(0.05, 0.04),
                 switch.prop = list(0.3, 0)),
    "h12.switch.hazard")
  # a per-cell list must have one element per subgroup cell
  expect_error(
    simdata_fast(nsim = 5, n = 50, a.time = c(0, 1), a.prop = 1,
                 h01.hazard = list(0.1, 0.2, 0.3), h02.hazard = 0.05,
                 prevalence = c(0.5, 0.5)),
    "per-cell list")
  # switching probabilities are checked with subgroups
  expect_error(
    simdata_fast(nsim = 5, n = 50, a.time = c(0, 1), a.prop = 1,
                 h01.hazard = 0.1, h02.hazard = 0.05, switch.prop = 1.5,
                 h12.switch.hazard = 0.02, prevalence = c(0.5, 0.5)),
    "probability")
  # only clock-reset is implemented
  expect_error(
    simdata_fast(nsim = 5, n = 50, a.time = c(0, 1), a.prop = 1,
                 h01.hazard = 0.1, h02.hazard = 0.05, switch.clock = "forward"),
    "reset")
  # the illness-death model and the single-endpoint argument are exclusive
  expect_error(
    simdata_fast(nsim = 5, n = 50, a.time = c(0, 1), a.prop = 1,
                 h01.hazard = 0.1, h02.hazard = 0.05, e.hazard = 0.1),
    "h01")
})

test_that("simdata_fast: single-endpoint path is not diverted by the new feature", {
  # Without any h01 / h02 argument the original single-endpoint schema is used.
  dat <- simdata_fast(
    nsim = 5, n = c(40, 40), a.time = c(0, 12), a.prop = 1,
    e.hazard = list(log(2) / 12, log(2) / 18), seed = 9
  )
  expect_true("surv_time" %in% names(dat))
  expect_false("e1_surv_time" %in% names(dat))
  expect_false("intermediate" %in% names(dat))
})

test_that("simdata_fast (illness-death): a per-group dropout list gives two groups", {
  dat <- simdata_fast(nsim = 2, n = 100, a.time = c(0, 10), a.rate = 10,
                      h01.hazard = 0.1, h02.hazard = 0.05,
                      d.hazard = list(0.01, 0.03), seed = 1)
  expect_equal(sort(unique(dat$group)), c(1, 2))
  expect_equal(nrow(dat), 200)
})

test_that("simdata_fast (illness-death): scalar n with alloc gives two groups", {
  dat <- simdata_fast(nsim = 2, n = 100, alloc = c(1, 1), a.time = c(0, 10),
                      a.rate = 10, h01.hazard = 0.1, h02.hazard = 0.05,
                      seed = 1)
  expect_equal(sort(unique(dat$group)), c(1, 2))
  expect_equal(sum(dat$group == 1 & dat$sim == 1), 50)
})

test_that("simdata_fast (illness-death): a switch requires an observed intermediate event", {
  dat <- simdata_fast(
    nsim = 1, n = c(20000, 20000), a.time = c(0, 1), a.prop = 1,
    h01.hazard = list(0.10, 0.07), h02.hazard = list(0.05, 0.04),
    d.hazard = 0.08, switch.prop = list(0.5, 0),
    h12.switch.hazard = list(0.015, 0.015), seed = 9
  )
  sw <- dat$switched == 1L
  expect_true(any(sw))
  expect_true(all(dat$intermediate[sw] == 1L & dat$e1_event[sw] == 1L))
  # A latent intermediate event after dropout does not lead to a switch.
  late <- dat$intermediate == 1L & dat$e1_event == 0L
  expect_true(any(late))
  expect_true(all(dat$switched[late] == 0L))
  expect_true(all(is.na(dat$switch_time[late])))
  obs <- dat$group == 1L & dat$intermediate == 1L & dat$e1_event == 1L
  expect_lt(abs(mean(dat$switched[obs]) - 0.5), 0.02)
})

test_that("simdata_fast (illness-death): infinite latent times are censored", {
  dat <- simdata_fast(
    nsim = 1, n = c(3000, 3000), a.time = c(0, 1), a.prop = 1,
    h01.hazard = c(0.05, 0), h01.time = c(0, 12, Inf),
    h02.hazard = c(0.02, 0), h02.time = c(0, 12, Inf), seed = 4
  )
  never <- !is.finite(dat$e1_surv_time)
  expect_true(any(never))
  expect_equal(dat$e1_event, as.integer(is.finite(dat$e1_surv_time)))
  expect_equal(dat$e2_event, as.integer(is.finite(dat$e2_surv_time)))
  expect_true(all(dat$intermediate[never] == 0L))
})

test_that("simdata_fast (illness-death): a single cell reproduces the data without subgroups", {
  args <- list(nsim = 5, n = c(60, 60), a.time = c(0, 6), a.rate = 20,
               h01.hazard = list(0.10, 0.07), h02.hazard = list(0.03, 0.02),
               h12.hazard = list(0.08, 0.08), switch.prop = list(0.4, 0),
               h12.switch.hazard = list(0.04, 0.04), d.hazard = 0.01,
               seed = 12)
  a <- do.call(simdata_fast, c(args, list(prevalence = 1)))
  b <- do.call(simdata_fast, args)
  expect_true(all(a$subgroup == 1L))
  expect_identical(a[names(b)], b)
  a1 <- simdata_fast(nsim = 3, n = 50, a.time = c(0, 5), a.rate = 10,
                     h01.hazard = 0.1, h02.hazard = 0.05, prevalence = 1,
                     seed = 3)
  b1 <- simdata_fast(nsim = 3, n = 50, a.time = c(0, 5), a.rate = 10,
                     h01.hazard = 0.1, h02.hazard = 0.05, seed = 3)
  expect_identical(a1[names(b1)], b1)
})

test_that("simdata_fast (illness-death): subgroups follow their own transition hazards", {
  h01 <- list(list(0.10, 0.05), list(0.07, 0.04))
  h02 <- c(0.03, 0.02)
  h12 <- list(c(0.08, 0.20), c(0.06, 0.06))
  dat <- simdata_fast(nsim = 1, n = c(40000, 40000), a.time = c(0, 1),
                      a.prop = 1, h01.hazard = h01,
                      h02.hazard = list(h02[1], h02[2]),
                      h12.hazard = list(list(0.08, 0.20), 0.06),
                      prevalence = c(0.3, 0.7), seed = 21)
  expect_lt(abs(mean(dat$subgroup == 1) - 0.3), 0.01)
  for (g in 1:2) {
    for (s in 1:2) {
      sel <- dat$group == g & dat$subgroup == s
      l01 <- h01[[g]][[s]]
      l02 <- h02[g]
      # Fleischer model: PFS is exponential with rate h01 + h02, the
      # intermediate event occurs with probability h01 / (h01 + h02), and the
      # post-event survival is exponential with rate h12 (clock-reset).
      expect_equal(1 / mean(dat$e1_surv_time[sel]), l01 + l02,
                   tolerance = 0.03)
      expect_lt(abs(mean(dat$intermediate[sel]) - l01 / (l01 + l02)), 0.02)
      prog <- sel & dat$intermediate == 1L
      post <- dat$e2_surv_time[prog] - dat$e1_surv_time[prog]
      expect_equal(1 / mean(post), h12[[g]][s], tolerance = 0.04)
    }
  }
})

test_that("simdata_fast (illness-death): a subgroup matches its separate simulation with scaled accrual", {
  # One group of 60,000 with 25 percent in subgroup 1, against 15,000
  # subjects of subgroup 1 simulated alone with the accrual rate scaled by the
  # prevalence. The distributions of accrual and of both endpoints agree.
  mix <- simdata_fast(nsim = 1, n = 60000, a.time = c(0, 12), a.rate = 5000,
                      h01.hazard = list(0.10, 0.04), h02.hazard = 0.03,
                      prevalence = c(0.25, 0.75), fixed.alloc = TRUE,
                      seed = 31)
  sep <- simdata_fast(nsim = 1, n = 15000, a.time = c(0, 12), a.rate = 1250,
                      h01.hazard = 0.10, h02.hazard = 0.03, seed = 32)
  s1 <- mix[mix$subgroup == 1L, ]
  expect_equal(nrow(s1), 15000L)
  expect_lt(abs(mean(s1$accrual_time) - mean(sep$accrual_time)), 0.15)
  for (t in c(3, 6, 9)) {
    expect_lt(abs(mean(s1$accrual_time <= t) - mean(sep$accrual_time <= t)),
              0.02)
  }
  for (t in c(6, 12, 24)) {
    expect_lt(abs(mean(s1$e1_surv_time > t) - mean(sep$e1_surv_time > t)),
              0.02)
    expect_lt(abs(mean(s1$e2_surv_time > t) - mean(sep$e2_surv_time > t)),
              0.02)
    expect_lt(abs(mean(s1$e2_calendar_time <= t + 6) -
                    mean(sep$e2_calendar_time <= t + 6)), 0.02)
  }
})

test_that("simdata_fast (illness-death): subgroup data work with the analysis functions", {
  dat <- simdata_fast(
    nsim = 20, n = c(150, 150), a.time = c(0, 12), a.rate = 25,
    h01.hazard = list(list(0.10, 0.08), list(0.07, 0.04)),
    h02.hazard = list(0.03, 0.02), d.hazard = 0.01,
    prevalence = list(control = c(0.5, 0.5), treatment = c(0.4, 0.6)),
    seed = 41
  )
  expect_identical(names(dat)[1:4], c("sim", "group", "subgroup",
                                      "accrual_time"))
  expect_lt(abs(mean(dat$subgroup[dat$group == 2] == 1) - 0.4), 0.05)
  pfs <- data.frame(sim = dat$sim, group = dat$group,
                    subgroup = dat$subgroup, accrual_time = dat$accrual_time,
                    tte = dat$e1_tte, event = dat$e1_event)
  res <- analysis_fast(pfs, control = 1, event.looks = 120, by.subgroup = TRUE)
  expect_setequal(unique(res$population),
                  c("overall", "subgroup_1", "subgroup_2"))
  cut <- cutoff_fast(dat, event.looks = 40, event.subset = dat$subgroup == 2,
                     tte.col = "e1_tte", event.col = "e1_event")
  expect_true(all(is.finite(cut[, 1])))
  sw <- switch_fast(dat, group = 1, prob = 0.5, when = "intermediate",
                    aft.factor = 1.3, seed = 1)
  expect_identical(sw$subgroup, dat$subgroup)
  expect_true(any(sw$switched == 1))
})

test_that("simdata_fast (illness-death): switching and dropout can differ by subgroup", {
  dat <- simdata_fast(
    nsim = 1, n = c(30000, 30000), a.time = c(0, 1), a.prop = 1,
    h01.hazard = list(0.10, 0.07), h02.hazard = list(0.03, 0.02),
    switch.prop = list(list(0.6, 0.2), 0),
    h12.switch.hazard = list(list(0.02, 0.05), 0.05),
    d.hazard = list(list(0.01, 0.05), 0.01),
    prevalence = c(0.5, 0.5), seed = 51
  )
  for (s in 1:2) {
    ctl <- dat$group == 1L & dat$subgroup == s
    obs <- ctl & dat$intermediate == 1L & dat$e1_event == 1L
    expect_lt(abs(mean(dat$switched[obs]) - c(0.6, 0.2)[s]), 0.02)
    sw <- ctl & dat$switched == 1L
    post <- dat$e2_surv_time[sw] - dat$e1_surv_time[sw]
    expect_equal(1 / mean(post), c(0.02, 0.05)[s], tolerance = 0.05)
    expect_equal(1 / mean(dat$dropout_time[ctl]), c(0.01, 0.05)[s],
                 tolerance = 0.03)
  }
  expect_true(all(dat$switched[dat$group == 2L] == 0L))
})

test_that("simdata_fast (illness-death): switch.prop and per-group lists are checked without subgroups", {
  # A percentage instead of a probability.
  expect_error(
    simdata_fast(nsim = 1, n = c(200, 200), a.time = c(0, 6), a.prop = 1,
                 h01.hazard = list(0.10, 0.07), h02.hazard = list(0.03, 0.02),
                 switch.prop = list(40, 0), h12.switch.hazard = 0.05, seed = 1),
    "single probability in \\[0, 1\\]")
  expect_error(
    simdata_fast(nsim = 1, n = c(200, 200), a.time = c(0, 6), a.prop = 1,
                 h01.hazard = list(0.10, 0.07), h02.hazard = list(0.03, 0.02),
                 switch.prop = list(-0.2, 0), h12.switch.hazard = 0.05,
                 seed = 1),
    "single probability in \\[0, 1\\]")
  # A per-group list with a third element.
  expect_error(
    simdata_fast(nsim = 1, n = c(50, 50), a.time = c(0, 6), a.prop = 1,
                 h01.hazard = list(0.10, 0.07, 0.05),
                 h02.hazard = list(0.03, 0.02), seed = 1),
    "list of length 2")
  # Valid probabilities are accepted.
  d <- simdata_fast(nsim = 1, n = c(50, 50), a.time = c(0, 6), a.prop = 1,
                    h01.hazard = list(0.10, 0.07), h02.hazard = list(0.03, 0.02),
                    switch.prop = list(0.4, 0), h12.switch.hazard = 0.05,
                    seed = 1)
  expect_equal(nrow(d), 100L)
})

test_that("simdata_fast (illness-death): n and nsim within rounding error are rounded", {
  ns <- 100 * 0.29
  n <- 90 * c(0.3, 0.7)  # 63 is 62.99999999999999
  d1 <- simdata_fast(nsim = ns, n = n, a.time = c(0, 12), a.prop = 1,
                     h01.hazard = list(0.10, 0.07),
                     h02.hazard = list(0.03, 0.02), seed = 1)
  tab1 <- table(d1$sim, d1$group)
  expect_equal(nrow(tab1), 29L)
  expect_true(all(tab1[, "1"] == 27L) && all(tab1[, "2"] == 63L))
  # With subgroups and fixed cell counts (9 and 18, 21 and 42 per trial).
  d2 <- simdata_fast(nsim = ns, n = n, a.time = c(0, 12), a.prop = 1,
                     h01.hazard = list(0.10, 0.07),
                     h02.hazard = list(0.03, 0.02),
                     prevalence = c(1, 2) / 3, fixed.alloc = TRUE, seed = 1)
  tab2 <- table(d2$group, d2$subgroup, d2$sim)
  expect_equal(dim(tab2)[3L], 29L)
  expect_equal(as.vector(tab2), rep(c(9, 21, 18, 42), 29L))
})

test_that("simdata_fast (illness-death): the list checks name the median arguments", {
  expect_error(
    simdata_fast(nsim = 1, n = c(50, 50), a.time = c(0, 6), a.prop = 1,
                 h01.median = list(7, 10, 14), h02.median = list(20, 30),
                 seed = 1),
    "'h01.hazard' (or 'h01.median')", fixed = TRUE)
  expect_error(
    simdata_fast(nsim = 1, n = c(50, 50), a.time = c(0, 6), a.prop = 1,
                 h01.hazard = list(0.10, 0.07), h02.hazard = list(0.03, 0.02),
                 d.median = list(60, 60, 60), seed = 1),
    "'d.hazard' (or 'd.median')", fixed = TRUE)
  # With subgroups: the per-group and the per-cell checks.
  expect_error(
    simdata_fast(nsim = 1, n = c(50, 50), a.time = c(0, 6), a.prop = 1,
                 h01.median = list(7, 10, 14), h02.median = list(20, 30),
                 prevalence = c(0.5, 0.5), seed = 1),
    "'h01.hazard' (or 'h01.median')", fixed = TRUE)
  expect_error(
    simdata_fast(nsim = 1, n = c(50, 50), a.time = c(0, 6), a.prop = 1,
                 h01.median = list(list(7, 10, 14), list(8, 12, 16)),
                 h02.median = list(20, 30), prevalence = c(0.5, 0.5),
                 seed = 1),
    "A per-cell list for 'h01.hazard' (or 'h01.median')", fixed = TRUE)
})
