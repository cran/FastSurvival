test_that("medsurv_fast returns the expected structure for two groups", {
  set.seed(1)
  n <- 200
  g <- rep(0:1, each = n / 2)
  tt <- rexp(n, rate = ifelse(g == 0, 0.1, 0.07))
  cc <- rexp(n, rate = 0.02)
  time <- pmin(tt, cc)
  event <- as.integer(tt <= cc)

  res <- medsurv_fast(time, event, group = g, control = 0)

  expect_s3_class(res, "medsurv_fast")
  expect_true(all(c("median.control", "median.treatment", "diff",
                    "se.diff", "z", "chisq", "p") %in% names(res)))
  expect_equal(unname(res["chisq"]), unname(res["z"])^2, tolerance = 1e-8)
  expect_true(is.finite(res["p"]) && res["p"] >= 0 && res["p"] <= 1)
})

test_that("median estimate matches survival::survfit", {
  skip_if_not_installed("survival")
  set.seed(2)
  n <- 300
  tt <- rexp(n, 0.08)
  cc <- rexp(n, 0.02)
  time <- pmin(tt, cc)
  event <- as.integer(tt <= cc)
  # Only the reference computation is protected; an error in medsurv_fast
  # fails the test instead of skipping it.
  est <- unname(medsurv_fast(time, event)["median"])
  ref <- tryCatch({
    sf <- survival::survfit(survival::Surv(time, event) ~ 1)
    unname(summary(sf)$table["median"])
  }, error = function(e) NULL)
  skip_if(is.null(ref) || !is.finite(ref), "survival comparison unavailable")
  expect_equal(est, ref, tolerance = 1e-6)
})

test_that("difference confidence interval matches diff plus or minus z times se", {
  set.seed(3)
  n <- 240
  g <- rep(0:1, each = n / 2)
  tt <- rexp(n, rate = ifelse(g == 0, 0.1, 0.06))
  cc <- rexp(n, rate = 0.015)
  time <- pmin(tt, cc)
  event <- as.integer(tt <= cc)

  res <- medsurv_fast(time, event, group = g, control = 0, conf.level = 0.95)
  zc <- qnorm(0.975)
  expect_equal(unname(res["lower.diff"]),
               unname(res["diff"] - zc * res["se.diff"]), tolerance = 1e-8)
  expect_equal(unname(res["upper.diff"]),
               unname(res["diff"] + zc * res["se.diff"]), tolerance = 1e-8)
})

test_that("one-sided and two-sided p-values are consistent with z", {
  set.seed(4)
  n <- 220
  g <- rep(0:1, each = n / 2)
  tt <- rexp(n, rate = ifelse(g == 0, 0.12, 0.06))
  cc <- rexp(n, rate = 0.015)
  time <- pmin(tt, cc)
  event <- as.integer(tt <= cc)

  res2 <- medsurv_fast(time, event, group = g, control = 0, side = 2)
  res1 <- medsurv_fast(time, event, group = g, control = 0, side = 1)
  z <- unname(res2["z"])
  expect_equal(unname(res1["p"]), pnorm(z, lower.tail = FALSE), tolerance = 1e-8)
  expect_equal(unname(res2["p"]), 2 * pnorm(-abs(z)), tolerance = 1e-8)
})

test_that("single-group mode returns a median and confidence interval without a test", {
  set.seed(5)
  n <- 200
  tt <- rexp(n, 0.1)
  cc <- rexp(n, 0.02)
  time <- pmin(tt, cc)
  event <- as.integer(tt <= cc)

  res <- medsurv_fast(time, event)
  expect_s3_class(res, "medsurv_fast")
  expect_true(all(c("median", "se", "lower", "upper") %in% names(res)))
  expect_false("p" %in% names(res))
  expect_true(unname(res["lower"]) <= unname(res["median"]))
  expect_true(unname(res["upper"]) >= unname(res["median"]))
})

test_that("bandwidth affects the standard error but not the median (km method)", {
  set.seed(6)
  n <- 240
  g <- rep(0:1, each = n / 2)
  tt <- rexp(n, rate = ifelse(g == 0, 0.1, 0.06))
  cc <- rexp(n, rate = 0.015)
  time <- pmin(tt, cc)
  event <- as.integer(tt <= cc)

  res_a <- medsurv_fast(time, event, group = g, control = 0, bw = 4)
  res_b <- medsurv_fast(time, event, group = g, control = 0, bw = 12)
  expect_equal(unname(res_a["median.control"]), unname(res_b["median.control"]))
  expect_equal(unname(res_a["median.treatment"]), unname(res_b["median.treatment"]))
  expect_false(isTRUE(all.equal(unname(res_a["se.control"]),
                                unname(res_b["se.control"]))))
})

test_that("method changes the standard error but not the median", {
  set.seed(202)
  n_per <- 250
  grp <- rep(0:1, each = n_per)
  tt <- c(rexp(n_per, log(2) / 12), rexp(n_per, log(2) / 16))
  cc <- rexp(2 * n_per, 0.01)
  time <- pmin(tt, cc)
  event <- as.integer(tt <= cc)

  km <- medsurv_fast(time, event, group = grp, control = 0, method = "km")
  np <- medsurv_fast(time, event, group = grp, control = 0, method = "nph")
  expect_equal(unname(km["median.control"]), unname(np["median.control"]))
  expect_equal(unname(km["median.treatment"]), unname(np["median.treatment"]))
  expect_false(isTRUE(all.equal(unname(km["se.diff"]), unname(np["se.diff"]))))
})

test_that("method = 'nph' reproduces nph::nphparams median inference", {
  skip_if_not_installed("nph")
  set.seed(101)
  n_per <- 300
  grp <- rep(0:1, each = n_per)
  tt <- c(rexp(n_per, log(2) / 12), rexp(n_per, log(2) / 16))
  cc <- rexp(2 * n_per, 0.01)
  time <- pmin(tt, cc)
  event <- as.integer(tt <= cc)
  # Only the reference computation is protected; an error in medsurv_fast
  # fails the test instead of skipping it.
  fast <- medsurv_fast(time, event, group = grp, control = 0, side = 2,
                       method = "nph")
  np <- tryCatch(
    nph::nphparams(time = time, event = event, group = grp,
                   param_type = "Q", param_par = 0.5),
    error = function(e) NULL)
  skip_if(is.null(np), "nph::nphparams failed")
  # An exact comparison requires the Kaplan-Meier and Nelson-Aalen medians to
  # coincide, since medsurv_fast keeps the Kaplan-Meier point estimate.
  skip_if(!isTRUE(all.equal(unname(fast["diff"]),
                            as.numeric(np$tab$Estimate), tolerance = 1e-8)),
          "Kaplan-Meier and Nelson-Aalen medians differ for this dataset")
  expect_equal(unname(fast["se.diff"]), as.numeric(np$tab$SE),
               tolerance = 1e-6)
  expect_equal(unname(fast["p"]), as.numeric(np$tab$p_unadj),
               tolerance = 1e-6)
})

test_that("input validation works", {
  expect_error(medsurv_fast(1:5, c(0, 1, 0, 1)), "same length")
  expect_error(medsurv_fast(1:5, rep(2L, 5)), "coded as 0")
  expect_error(
    medsurv_fast(1:6, rep(0:1, 3), group = rep(1:3, 2), control = 1),
    "two distinct"
  )
  expect_error(
    medsurv_fast(1:6, rep(0:1, 3), group = rep(0:1, 3)),
    "control must be specified"
  )
  expect_error(
    medsurv_fast(1:6, rep(0:1, 3), group = rep(0:1, 3), control = 0,
                 method = "foo"),
    "should be one of"
  )
})

test_that("type I error is approximately controlled under the null", {
  set.seed(42)
  nsim <- 200L
  reject <- 0L
  for (s in seq_len(nsim)) {
    n <- 200
    g <- rep(0:1, each = n / 2)
    tt <- rexp(n, rate = 0.1)
    cc <- rexp(n, rate = 0.02)
    time <- pmin(tt, cc)
    event <- as.integer(tt <= cc)
    res <- medsurv_fast(time, event, group = g, control = 0, side = 2)
    if (is.finite(res["p"]) && res["p"] < 0.05) reject <- reject + 1L
  }
  expect_lt(reject / nsim, 0.12)
})

test_that("the median follows the survfit convention on a flat stretch at 0.5", {
  # Hand calculation: uncensored times 1..10 give S = 0.5 on [5, 6), so the
  # median is the midpoint 5.5.
  expect_equal(unname(medsurv_fast(1:10, rep(1L, 10))["median"]), 5.5)
  # Times 1..24: S(12) = 0.5 up to floating-point rounding, so the median is
  # 12.5 rather than 13.
  expect_equal(unname(medsurv_fast(1:24, rep(1L, 24))["median"]), 12.5)
  # Times 1..5: S drops from 0.6 to 0.4 at t = 3, so the median is 3.
  expect_equal(unname(medsurv_fast(1:5, rep(1L, 5))["median"]), 3)
})

test_that("medsurv_fast: presorted = TRUE checks the order", {
  time  <- c(3, 1, 2, 4, 6, 5)
  event <- c(1, 1, 0, 1, 1, 0)
  group <- c(0, 1, 0, 1, 0, 1)
  expect_error(medsurv_fast(time, event, group, control = 0, presorted = TRUE), "presorted = FALSE")
})

test_that("medsurv_fast: a median at 0.5 up to an infinite time is NA", {
  # In the treatment group 5 of 10 subjects have events at times 1 to 5, so the
  # Kaplan-Meier estimate is exactly 0.5 after time 5, and the other 5 have an
  # infinite time (no finite event or dropout time). The midpoint of the flat
  # stretch at 0.5 is not defined. The control group reaches 0.5 at time 5 and
  # its next event is at time 6, so its median is 5.5.
  time  <- c(1:10, 1:5, rep(Inf, 5))
  event <- c(rep(1, 10), rep(1, 5), rep(0, 5))
  group <- rep(0:1, each = 10)
  for (m in c("km", "nph")) {
    res <- medsurv_fast(time, event, group = group, control = 0, method = m)
    expect_equal(unname(res["median.control"]), 5.5)
    for (nm in c("median.treatment", "diff", "se.diff", "z", "p")) {
      expect_true(is.na(res[[nm]]), info = paste(m, nm))
    }
  }
  one <- medsurv_fast(time[group == 1], event[group == 1])
  expect_true(is.na(one[["median"]]))
  # The same data with the infinite times replaced by a finite censoring time
  # give a finite median (the midpoint of 5 and 20).
  time_f <- time
  time_f[is.infinite(time_f)] <- 20
  res_f <- medsurv_fast(time_f, event, group = group, control = 0)
  expect_equal(unname(res_f["median.treatment"]), 12.5)
})

test_that("medsurv_fast: time and event must be numeric and time non-negative", {
  set.seed(1)
  tt <- rexp(40, 0.1)
  ee <- rbinom(40, 1, 0.7)
  expect_error(medsurv_fast(as.character(tt), ee), "numeric")
  expect_error(medsurv_fast(c(-1, tt[-1]), ee), "non-negative")
  expect_error(medsurv_fast(tt, factor(ee)), "numeric or logical")
  # A logical event indicator is accepted and gives the same result.
  expect_equal(unclass(medsurv_fast(tt, ee == 1)), unclass(medsurv_fast(tt, ee)))
})

test_that("medsurv_fast: the formula interface gives the vector result", {
  set.seed(67)
  n <- 120
  d <- data.frame(tt = round(rexp(n, 0.1), 1), ee = rbinom(n, 1, 0.8),
                  g = rep(c(1, 0), n / 2))
  for (m in c("km", "nph")) {
    expect_identical(
      medsurv_fast(Surv(tt, ee) ~ g, data = d, method = m),
      medsurv_fast(d$tt, d$ee, group = d$g, control = 0, method = m))
  }
  ref1 <- medsurv_fast(d$tt, d$ee)
  expect_identical(medsurv_fast(Surv(tt, ee) ~ 1, data = d), ref1)
  expect_identical(class(ref1), c("medsurv_fast", "fastsurv"))
  expect_error(medsurv_fast(d$tt, d$ee, data = d), "only with a formula")
})

test_that("medsurv_fast: the formula interface passes every argument", {
  set.seed(60)
  n <- 120
  d <- data.frame(tt = round(rexp(n, 0.1), 1) + 0.1, ee = rbinom(n, 1, 0.8),
                  g = rep(c("b", "a"), n / 2), s = rep(c("x", "y", "z"), n / 3))
  d1 <- d[order(d$tt), ]
  d2 <- d[order(d$s, d$tt), ]
  expect_identical(
    medsurv_fast(Surv(tt, ee) ~ g, d1, control = "b", side = 1,
                 conf.level = 0.9, conf.type = "plain", method = "nph",
                 bw = c(2, 3), presorted = TRUE),
    medsurv_fast(d1$tt, d1$ee, d1$g, "b", side = 1, conf.level = 0.9,
                 conf.type = "plain", method = "nph", bw = c(2, 3),
                 presorted = TRUE))
})
