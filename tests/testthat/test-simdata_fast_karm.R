test_that("simdata_fast: multi-arm output has the expected structure", {
  nsim  <- 50
  n_arm <- c(40, 50, 60)
  df <- simdata_fast(
    nsim     = nsim,
    n        = n_arm,
    a.time   = c(0, 12),
    a.rate   = sum(n_arm) / 12,
    e.median = list(12, 16, 20),
    seed     = 101
  )

  # One block of sum(n_arm) subjects per simulation.
  expect_equal(nrow(df), nsim * sum(n_arm))
  expect_true(all(c("sim", "group", "accrual_time", "surv_time",
                    "dropout_time", "tte", "event", "calendar_time") %in%
                    names(df)))

  # Arm labels are 1..K in the order of 'n'.
  expect_equal(sort(unique(df$group)), 1:3)

  # Per-arm counts within each simulation match n_arm exactly.
  tab <- table(df$sim, df$group)
  expect_true(all(tab[, 1] == n_arm[1]))
  expect_true(all(tab[, 2] == n_arm[2]))
  expect_true(all(tab[, 3] == n_arm[3]))

  # Rows are grouped by simulation (non-decreasing sim), as the analysis path
  # relies on for its skip-sort fast path.
  expect_false(is.unsorted(df$sim))
})

test_that("simdata_fast: multi-arm is reproducible from the seed", {
  args <- list(
    nsim     = 30,
    n        = c(30, 30, 30),
    a.time   = c(0, 12),
    a.rate   = 90 / 12,
    e.median = list(12, 16, 20),
    seed     = 202
  )
  df1 <- do.call(simdata_fast, args)
  df2 <- do.call(simdata_fast, args)
  expect_identical(df1, df2)

  # A different seed changes the realized data.
  args3 <- args; args3$seed <- 303
  df3 <- do.call(simdata_fast, args3)
  expect_false(identical(df1$surv_time, df3$surv_time))
})

test_that("simdata_fast: multi-arm recovers the per-arm hazards", {
  n_arm <- c(150, 150, 150)
  med   <- c(12, 16, 20)
  df <- simdata_fast(
    nsim     = 200,
    n        = n_arm,
    a.time   = c(0, 12),
    a.rate   = sum(n_arm) / 12,
    e.median = as.list(med),
    seed     = 404
  )
  # Exponential mean of the raw survival time is median / log(2), so the
  # per-arm empirical hazard 1 / mean(surv_time) recovers log(2) / median.
  for (g in 1:3) {
    est <- 1 / mean(df$surv_time[df$group == g])
    expect_equal(est, log(2) / med[g], tolerance = 0.03,
                 info = paste("arm", g))
  }
})

test_that("simdata_fast: multi-arm accepts e.hazard and per-arm dropout", {
  n_arm <- c(40, 40, 40)
  df <- simdata_fast(
    nsim     = 40,
    n        = n_arm,
    a.time   = c(0, 12),
    a.rate   = sum(n_arm) / 12,
    e.hazard = list(0.10, 0.08, 0.06),
    d.hazard = list(0.01, 0.01, 0.02),
    seed     = 505
  )
  expect_equal(nrow(df), 40 * sum(n_arm))
  expect_equal(sort(unique(df$group)), 1:3)
  # Dropout was requested, so some non-event censorings occur.
  expect_true(any(df$dropout_time < df$surv_time))
})

test_that("simdata_fast: multi-arm output feeds the pairwise analysis path", {
  nsim  <- 100
  n_arm <- c(120, 120, 120)
  df <- simdata_fast(
    nsim     = nsim,
    n        = n_arm,
    a.time   = c(0, 12),
    a.rate   = sum(n_arm) / 12,
    e.median = list(12, 16, 20),
    seed     = 606
  )
  # Control (group 1) versus experimental arm 2.
  sub <- df[df$group %in% c(1, 2), ]
  res <- analysis_fast(sub, control = 1, time.looks = 30,
                       stat = "logrank", side = 1)
  expect_equal(nrow(res), nsim)
  expect_true(all(is.finite(res$logrank.z)))
  expect_true(all(res$logrank.p >= 0 & res$logrank.p <= 1))
})

test_that("simdata_fast: multi-arm rejects unsupported combinations", {
  base_args <- list(
    nsim     = 10,
    n        = c(20, 20, 20),
    a.time   = c(0, 12),
    a.rate   = 60 / 12,
    e.median = list(12, 16, 20),
    seed     = 707
  )

  # Subgroups are not supported in multi-arm mode.
  bad_prev <- base_args
  bad_prev$prevalence <- c(0.5, 0.5)
  expect_error(do.call(simdata_fast, bad_prev), "does not support 'prevalence'")

  # The illness-death model is two-group only.
  bad_id <- base_args
  bad_id$h01.median <- list(8, 10, 12)
  bad_id$h02.median <- list(24, 28, 32)
  expect_error(do.call(simdata_fast, bad_id))
})

test_that("simdata_fast: length(n) > 2 without a per-arm list keeps the old error", {
  # A scalar hazard with a length-3 'n' is not multi-arm mode; the historical
  # validation error is preserved, so existing behavior is unchanged.
  expect_error(
    simdata_fast(nsim = 5, n = c(1, 2, 3), a.time = c(0, 1), a.prop = 1,
                 e.hazard = log(2) / 12),
    "scalar .total N. or a vector of length 2")
})
