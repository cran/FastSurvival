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
  # subgroups are not supported
  expect_error(
    simdata_fast(nsim = 5, n = 50, a.time = c(0, 1), a.prop = 1,
                 h01.hazard = 0.1, h02.hazard = 0.05, prevalence = c(0.5, 0.5)),
    "prevalence")
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
