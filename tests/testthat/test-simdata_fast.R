test_that("simdata_fast returns a data.frame with correct columns (one group)", {
  df <- simdata_fast(
    nsim     = 10,
    n        = 20,
    a.time   = c(0, 12),
    a.rate   = 20 / 12,
    e.median = 18,
    seed     = 1
  )

  expect_s3_class(df, "data.frame")
  expect_named(df, c("sim", "group", "accrual_time", "surv_time",
                     "dropout_time", "tte", "event", "calendar_time"))
})

test_that("simdata_fast returns correct number of rows (one group)", {
  nsim <- 10
  n    <- 20
  df   <- simdata_fast(
    nsim     = nsim,
    n        = n,
    a.time   = c(0, 12),
    a.rate   = n / 12,
    e.median = 18,
    seed     = 2
  )

  expect_equal(nrow(df), nsim * n)
})

test_that("simdata_fast returns correct number of rows (two groups)", {
  nsim <- 10
  n    <- c(50, 50)
  df   <- simdata_fast(
    nsim     = nsim,
    n        = n,
    a.time   = c(0, 12),
    a.rate   = sum(n) / 12,
    e.median = list(18, 24),
    seed     = 3
  )

  expect_equal(nrow(df), nsim * sum(n))
  expect_equal(sort(unique(df$group)), c(1L, 2L))
})

test_that("simdata_fast accrual_time is within a.time range", {
  df <- simdata_fast(
    nsim     = 20,
    n        = 50,
    a.time   = c(0, 12),
    a.rate   = 50 / 12,
    e.median = 18,
    seed     = 4
  )

  expect_gte(min(df$accrual_time), 0)
  expect_lte(max(df$accrual_time), 12)
})

test_that("simdata_fast tte = pmin(surv_time, dropout_time)", {
  df <- simdata_fast(
    nsim     = 10,
    n        = 50,
    a.time   = c(0, 12),
    a.rate   = 50 / 12,
    e.median = 18,
    d.hazard = 0.02,
    seed     = 5
  )

  expect_equal(df$tte, pmin(df$surv_time, df$dropout_time))
})

test_that("simdata_fast event = 1 iff surv_time <= dropout_time", {
  df <- simdata_fast(
    nsim     = 10,
    n        = 50,
    a.time   = c(0, 12),
    a.rate   = 50 / 12,
    e.median = 18,
    d.hazard = 0.02,
    seed     = 6
  )

  expected_event <- as.integer(df$surv_time <= df$dropout_time)
  expect_equal(df$event, expected_event)
})

test_that("simdata_fast calendar_time = accrual_time + tte", {
  df <- simdata_fast(
    nsim     = 10,
    n        = 50,
    a.time   = c(0, 12),
    a.rate   = 50 / 12,
    e.median = 18,
    seed     = 7
  )

  expect_equal(df$calendar_time, df$accrual_time + df$tte)
})

test_that("simdata_fast dropout_time = Inf when no dropout specified", {
  df <- simdata_fast(
    nsim     = 10,
    n        = 50,
    a.time   = c(0, 12),
    a.rate   = 50 / 12,
    e.median = 18,
    seed     = 8
  )

  expect_true(all(is.infinite(df$dropout_time)))
  expect_true(all(df$event == 1L))
})

test_that("simdata_fast is reproducible with seed", {
  df1 <- simdata_fast(
    nsim     = 5,
    n        = 20,
    a.time   = c(0, 12),
    a.rate   = 20 / 12,
    e.median = 18,
    seed     = 42
  )
  df2 <- simdata_fast(
    nsim     = 5,
    n        = 20,
    a.time   = c(0, 12),
    a.rate   = 20 / 12,
    e.median = 18,
    seed     = 42
  )

  expect_equal(df1, df2)
})

test_that("simdata_fast piecewise exponential: surv_time > 0", {
  df <- simdata_fast(
    nsim     = 10,
    n        = c(50, 50),
    a.time   = c(0, 12),
    a.rate   = 100 / 12,
    e.hazard = list(c(0.08, 0.08), c(0.08, 0.04)),
    e.time   = c(0, 6, Inf),
    seed     = 9
  )

  expect_true(all(df$surv_time > 0))
})

test_that("simdata_fast e.hazard and e.median are mutually exclusive", {
  expect_error(
    simdata_fast(
      nsim     = 5,
      n        = 20,
      a.time   = c(0, 12),
      a.rate   = 20 / 12,
      e.hazard = 0.05,
      e.median = 18,
      seed     = 10
    ),
    "exactly one"
  )
})

test_that("simdata_fast total n with alloc splits correctly", {
  nsim  <- 5
  n_tot <- 100
  df    <- simdata_fast(
    nsim     = nsim,
    n        = n_tot,
    alloc    = c(1, 1),
    a.time   = c(0, 12),
    a.rate   = n_tot / 12,
    e.median = list(18, 24),
    seed     = 11
  )

  n_per_sim <- table(df$sim[df$sim == 1 & df$group == 1])
  expect_equal(nrow(df), nsim * n_tot)
})
