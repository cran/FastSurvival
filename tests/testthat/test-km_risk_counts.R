test_that("km_risk_counts counts the subjects at risk (hand calculation)", {
  # Hand counts of the subjects with an observed time at or after each time
  time  <- c(1, 2, 2, 3, 4, 5, 0.5, 2, 6, 7)
  event <- c(1, 1, 0, 1, 0, 1, 1, 0, 1, 1)
  group <- rep(c(0, 1), c(6, 4))
  fit <- kmcurve_fast(time, event, group, control = 0)
  cnt <- km_risk_counts(fit, c(0, 2, 2.5, 5, 6, 8))
  expect_identical(rownames(cnt), c("control", "treat"))
  expect_equal(unname(cnt["control", ]), c(6, 5, 3, 1, 0, 0))
  expect_equal(unname(cnt["treat", ]), c(4, 3, 2, 2, 2, 0))
  # The control group is the one named by 'control', not the first value
  fit_rev <- kmcurve_fast(time, event, group, control = 1)
  expect_equal(unname(km_risk_counts(fit_rev, c(0, 2))["control", ]), c(4, 3))
  one <- km_risk_counts(kmcurve_fast(c(1, 2, 2, 3), c(1, 1, 0, 1)),
                        c(0, 2, 3.5))
  expect_equal(dim(one), c(1L, 3L))
  expect_equal(unname(one[1L, ]), c(4, 3, 0))
})

test_that("km_risk_counts matches survival::survfit", {
  set.seed(97)
  n  <- 80
  t0 <- rexp(2 * n, log(2) / 12)
  cc <- runif(2 * n, 0, 30)
  time  <- round(pmin(t0, cc), 1)
  event <- as.integer(t0 <= cc)
  group <- rep(c(1, 2), each = n)
  times <- c(0, 5, 10, 12.3, 20)
  cnt <- km_risk_counts(kmcurve_fast(time, event, group, control = 1), times)
  # Report whether the external comparison below is run or skipped.
  message("survival available for the number-at-risk comparison: ",
          requireNamespace("survival", quietly = TRUE))
  skip_if_not_installed("survival")
  for (g in 1:2) {
    ref <- survival::survfit(survival::Surv(time[group == g],
                                            event[group == g]) ~ 1)
    nr <- summary(ref, times = times, extend = TRUE)$n.risk
    expect_equal(unname(cnt[g, ]), nr)
  }
})
