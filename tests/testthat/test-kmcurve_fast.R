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
  med <- curve_median(c(0, km$te), c(1, km$surv))

  ref <- tryCatch(
    survival::survfit(survival::Surv(time, event) ~ 1),
    error = function(e) NULL)
  skip_if(is.null(ref), "survfit failed")
  ref_med <- unname(summary(ref)$table["median"])
  skip_if(is.na(ref_med), "median not reached")

  expect_equal(med, ref_med, tolerance = 0.5)
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

  expect_error(print(fit2), NA)
  expect_error(print(fit1), NA)
  expect_error(plot(fit2), NA)
  expect_error(
    plot(fit2, hr = TRUE, rmst = TRUE, tau = 18,
         conf.type = "log-log", ahr_line = TRUE),
    NA)
  expect_error(plot(fit1, rmst = TRUE, tau = 18, conf.type = "plain"), NA)
  expect_warning(plot(fit1, hr = TRUE))   # hazard ratio ignored for one group
})
