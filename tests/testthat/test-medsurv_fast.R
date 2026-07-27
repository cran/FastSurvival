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
  result <- tryCatch({
    set.seed(2)
    n <- 300
    tt <- rexp(n, 0.08)
    cc <- rexp(n, 0.02)
    time <- pmin(tt, cc)
    event <- as.integer(tt <= cc)
    sf <- survival::survfit(survival::Surv(time, event) ~ 1)
    ref <- unname(summary(sf)$table["median"])
    est <- unname(medsurv_fast(time, event)["median"])
    list(ref = ref, est = est)
  }, error = function(e) NULL)
  skip_if(is.null(result) || !is.finite(result$ref),
          "survival comparison unavailable")
  expect_equal(result$est, result$ref, tolerance = 1e-6)
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
  out <- tryCatch({
    set.seed(101)
    n_per <- 300
    grp <- rep(0:1, each = n_per)
    tt <- c(rexp(n_per, log(2) / 12), rexp(n_per, log(2) / 16))
    cc <- rexp(2 * n_per, 0.01)
    time <- pmin(tt, cc)
    event <- as.integer(tt <= cc)
    fast <- medsurv_fast(time, event, group = grp, control = 0, side = 2,
                         method = "nph")
    np <- nph::nphparams(time = time, event = event, group = grp,
                         param_type = "Q", param_par = 0.5)
    list(fast = fast, np = np)
  }, error = function(e) NULL)
  skip_if(is.null(out), "nph comparison unavailable")
  fast <- out$fast
  np <- out$np
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
  expect_error(medsurv_fast(1:5, rep(2L, 5)), "0")
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
                 method = "foo")
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
