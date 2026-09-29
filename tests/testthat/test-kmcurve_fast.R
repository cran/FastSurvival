test_that("kmcurve_fast builds a correct two-group object", {
  set.seed(11)
  n  <- 80
  t0 <- rexp(n, log(2) / 12)
  t1 <- rexp(n, log(2) / 18)
  cc <- runif(2 * n, 0, 30)
  time  <- pmin(c(t0, t1), cc)
  event <- as.integer(c(t0, t1) <= cc)
  group <- rep(c(1, 2), each = n)

  fit <- kmcurve_fast(time, event, group, control = 1)

  expect_s3_class(fit, "kmcurve_fast")
  expect_true(fit$two_group)
  expect_equal(unname(fit$n[["control"]]), n)
  expect_equal(unname(fit$n[["treat"]]), n)
  expect_equal(unname(fit$events[["control"]]), sum(event[group == 1]))
  expect_equal(unname(fit$events[["treat"]]),   sum(event[group == 2]))

  s <- fit$km$control$surv
  expect_true(all(s >= 0 & s <= 1))
  expect_true(all(diff(s) <= 1e-12))
})

test_that("kmcurve_fast supports the single-group case", {
  set.seed(12)
  n  <- 100
  t0 <- rexp(n, log(2) / 12)
  cc <- runif(n, 0, 30)
  time  <- pmin(t0, cc)
  event <- as.integer(t0 <= cc)

  f_null <- kmcurve_fast(time, event)                 # group = NULL
  f_one  <- kmcurve_fast(time, event, rep(1, n))      # single distinct value

  expect_false(f_null$two_group)
  expect_false(f_one$two_group)
  expect_null(f_null$km$treat)
  expect_equal(unname(f_null$n[["control"]]), n)
  expect_equal(f_null$km$control$surv, f_one$km$control$surv)
})

test_that("kmcurve_fast validates its arguments", {
  time  <- 1:6
  event <- c(1, 0, 1, 1, 0, 1)
  g3    <- c(1, 1, 2, 2, 3, 3)
  g2    <- c(1, 1, 1, 2, 2, 2)

  expect_error(kmcurve_fast(time, event, g3), "at most two")
  expect_error(kmcurve_fast(time, event, g2), "control")
  expect_error(kmcurve_fast(time, event, g2, control = 9), "control")
  expect_error(kmcurve_fast(time, event[1:3], g2), "same length")
})

test_that("kmcurve_fast matches survival::survfit (survival and SE)", {
  skip_if_not_installed("survival")
  set.seed(21)
  n  <- 120
  t0 <- rexp(n, log(2) / 12)
  t1 <- rexp(n, log(2) / 20)
  cc <- runif(2 * n, 0, 36)
  time  <- pmin(c(t0, t1), cc)
  event <- as.integer(c(t0, t1) <= cc)
  group <- rep(c(1, 2), each = n)

  fit <- kmcurve_fast(time, event, group, control = 1)

  ref <- tryCatch(
    survival::survfit(
      survival::Surv(time[group == 1], event[group == 1]) ~ 1,
      conf.type = "log"),
    error = function(e) NULL)
  skip_if(is.null(ref), "survfit failed")

  te <- fit$km$control$te
  sm <- summary(ref, times = te)

  expect_equal(fit$km$control$surv, sm$surv,    tolerance = 1e-8)
  expect_equal(fit$km$control$se,   sm$std.err, tolerance = 1e-8)
})

test_that("km_ci matches survival::survfit for log and log-log", {
  skip_if_not_installed("survival")
  set.seed(22)
  n  <- 120
  t0 <- rexp(n, log(2) / 12)
  cc <- runif(n, 0, 36)
  time  <- pmin(t0, cc)
  event <- as.integer(t0 <= cc)

  fit <- kmcurve_fast(time, event)
  km  <- fit$km$control
  te  <- km$te
  keep <- seq_len(length(te) - 1L)   # drop degenerate tail point

  for (ct in c("log", "log-log")) {
    ref <- tryCatch(
      survival::survfit(survival::Surv(time, event) ~ 1, conf.type = ct),
      error = function(e) NULL)
    skip_if(is.null(ref), "survfit failed")

    sm <- summary(ref, times = te)
    ci <- km_ci(km$surv, km$se, 0.95, ct)

    expect_equal(ci$lower[keep], sm$lower[keep], tolerance = 1e-6)
    expect_equal(ci$upper[keep], sm$upper[keep], tolerance = 1e-6)
  }
})

test_that("median read from the curve matches survival::survfit", {
  skip_if_not_installed("survival")
  set.seed(23)
  n  <- 200
  t0 <- rexp(n, log(2) / 12)
  cc <- runif(n, 0, 48)
  time  <- pmin(t0, cc)
  event <- as.integer(t0 <= cc)

  fit <- kmcurve_fast(time, event)
  km  <- fit$km$control
  med <- km_step_median(km$te, km$surv, km$tmax)

  ref <- tryCatch(
    survival::survfit(survival::Surv(time, event) ~ 1),
    error = function(e) NULL)
  skip_if(is.null(ref), "survfit failed")
  ref_med <- unname(summary(ref)$table["median"])
  skip_if(is.na(ref_med), "median not reached")

  expect_equal(med, ref_med, tolerance = 1e-10)
})

test_that("the printed median follows the survfit step convention", {
  # Uncensored times 1..10: S(5) = 0.5 exactly on [5, 6), so survfit reports
  # the midpoint 5.5 (hand calculation).
  fit <- kmcurve_fast(1:10, rep(1L, 10))
  km  <- fit$km$control
  expect_equal(km_step_median(km$te, km$surv, km$tmax), 5.5)
  # Uncensored times 1..24: the product of the Kaplan-Meier factors at t = 12
  # is 0.5 up to rounding error, so the midpoint 12.5 is expected.
  fit24 <- kmcurve_fast(1:24, rep(1L, 24))
  km24  <- fit24$km$control
  expect_equal(km_step_median(km24$te, km24$surv, km24$tmax), 12.5)
  # S drops from 0.6 at t = 2 to 0.4 at t = 3 (n = 5, no ties): the median is
  # the event time 3, not an interpolated value.
  fit5 <- kmcurve_fast(1:5, rep(1L, 5))
  km5  <- fit5$km$control
  expect_equal(km_step_median(km5$te, km5$surv, km5$tmax), 3)
})

test_that("plot and print methods run for both single and two-group objects", {
  set.seed(31)
  n  <- 60
  t0 <- rexp(n, log(2) / 12)
  t1 <- rexp(n, log(2) / 18)
  cc <- runif(2 * n, 0, 30)
  time  <- pmin(c(t0, t1), cc)
  event <- as.integer(c(t0, t1) <= cc)
  group <- rep(c(1, 2), each = n)

  fit2 <- kmcurve_fast(time, event, group, control = 1)
  fit1 <- kmcurve_fast(time[group == 1], event[group == 1])

  grDevices::pdf(tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)

  expect_output(print(fit2), "two groups")
  expect_output(print(fit1), "one group")
  expect_error(plot(fit2), NA)
  expect_error(
    plot(fit2, hr = TRUE, rmst = TRUE, tau = 18,
         conf.type = "log-log", ahr_line = TRUE),
    NA)
  expect_error(plot(fit1, rmst = TRUE, tau = 18, conf.type = "plain"), NA)
  expect_warning(plot(fit1, hr = TRUE))   # hazard ratio ignored for one group
})
