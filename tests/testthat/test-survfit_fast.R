test_that("survfit_fast returns named vector of length 4", {
  set.seed(1)
  t_raw <- rexp(50, rate = 0.1)
  e_raw <- rbinom(50, 1, 0.7)
  ord   <- order(t_raw)
  res   <- survfit_fast(t_raw[ord], e_raw[ord], t_eval = 8)

  expect_length(res, 4L)
  expect_named(res, c("surv", "std.err", "lower", "upper"))
})

test_that("survfit_fast surv is in [0, 1]", {
  set.seed(2)
  t_raw <- rexp(100, rate = 0.1)
  e_raw <- rbinom(100, 1, 0.7)
  ord   <- order(t_raw)

  for (t_eval in c(1, 5, 10, 20)) {
    res <- survfit_fast(t_raw[ord], e_raw[ord], t_eval = t_eval)
    expect_gte(unname(res["surv"]), 0)
    expect_lte(unname(res["surv"]), 1)
  }
})

test_that("survfit_fast CI lower <= surv <= upper", {
  set.seed(3)
  t_raw <- rexp(100, rate = 0.1)
  e_raw <- rbinom(100, 1, 0.8)
  ord   <- order(t_raw)

  for (type in c("plain", "log", "log-log")) {
    res <- survfit_fast(t_raw[ord], e_raw[ord], t_eval = 8, conf.type = type)
    if (!any(is.na(res))) {
      expect_lte(unname(res["lower"]), unname(res["surv"]) + 1e-10)
      expect_gte(unname(res["upper"]), unname(res["surv"]) - 1e-10)
    }
  }
})

test_that("survfit_fast returns surv=1 when t_eval before first observation", {
  set.seed(4)
  t_raw <- rexp(50, rate = 0.1) + 5  # all times > 5
  e_raw <- rep(1L, 50)
  ord   <- order(t_raw)
  res   <- survfit_fast(t_raw[ord], e_raw[ord], t_eval = 0.01)

  expect_equal(unname(res["surv"]),    1)
  expect_equal(unname(res["std.err"]), 0)
  expect_equal(unname(res["lower"]),   1)
  expect_equal(unname(res["upper"]),   1)
})

test_that("survfit_fast returns all-NA when n=0", {
  res <- survfit_fast(numeric(0), integer(0), t_eval = 5)
  expect_true(all(is.na(res)))
})

test_that("survfit_fast presorted=TRUE and presorted=FALSE agree", {
  set.seed(5)
  t_raw <- rexp(100, rate = 0.1)
  e_raw <- rbinom(100, 1, 0.7)
  ord   <- order(t_raw)

  res_pre <- survfit_fast(t_raw[ord], e_raw[ord], t_eval = 10,
                          presorted = TRUE)
  res_uns <- survfit_fast(t_raw, e_raw, t_eval = 10,
                          presorted = FALSE)
  expect_equal(unname(res_pre), unname(res_uns), tolerance = 1e-12)
})

test_that("survfit_fast agrees with survival::survfit (plain CI)", {
  skip_if_not_installed("survival")
  set.seed(6)
  t_raw  <- rexp(100, rate = 0.1)
  e_raw  <- rbinom(100, 1, 0.7)
  t_eval <- 8

  ord <- order(t_raw)
  res <- survfit_fast(t_raw[ord], e_raw[ord], t_eval = t_eval,
                      conf.type = "plain")

  fit <- survival::survfit(survival::Surv(t_raw, e_raw) ~ 1,
                           conf.type = "plain")
  ref <- summary(fit, times = t_eval, extend = TRUE)

  expect_equal(unname(res["surv"]),  ref$surv,  tolerance = 1e-6)
  expect_equal(unname(res["lower"]), ref$lower, tolerance = 1e-6)
  expect_equal(unname(res["upper"]), ref$upper, tolerance = 1e-6)
})

test_that("survfit_fast agrees with survival::survfit (log CI)", {
  skip_if_not_installed("survival")
  set.seed(7)
  t_raw  <- rexp(100, rate = 0.1)
  e_raw  <- rbinom(100, 1, 0.7)
  t_eval <- 10

  ord <- order(t_raw)
  res <- survfit_fast(t_raw[ord], e_raw[ord], t_eval = t_eval,
                      conf.type = "log")

  fit <- survival::survfit(survival::Surv(t_raw, e_raw) ~ 1,
                           conf.type = "log")
  ref <- summary(fit, times = t_eval, extend = TRUE)

  expect_equal(unname(res["surv"]),  ref$surv,  tolerance = 1e-6)
  expect_equal(unname(res["lower"]), ref$lower, tolerance = 1e-6)
  expect_equal(unname(res["upper"]), ref$upper, tolerance = 1e-6)
})

test_that("survfit_fast caps the log interval at 1 and validates its input", {
  # Hand calculation: S(1) = 0.99 with SE 0.99 * sqrt(1 / (100 * 99)), so the
  # unrestricted upper limit 0.99 * exp(1.96 * 0.01005) is about 1.0097.
  res <- survfit_fast(1:100, rep(1L, 100), t_eval = 1, conf.type = "log")
  expect_equal(unname(res["surv"]), 0.99)
  expect_equal(unname(res["upper"]), 1)
  expect_error(survfit_fast(c(3, 1, 2), c(1, 1, 1), t_eval = 2), "not sorted")
  expect_error(survfit_fast(1:3, c(1, 1), t_eval = 2), "same length")
  expect_error(survfit_fast(1:3, c(1, 2, 1), t_eval = 2), "coded as 0")
  expect_error(survfit_fast(1:3, c(1, 1, 1), t_eval = 2, conf.level = 95),
               "conf.level")
})

test_that("survfit_fast: negative times are rejected", {
  expect_error(survfit_fast(c(-1, 2, 3), c(1, 0, 1), t_eval = 2,
                            presorted = FALSE), "non-negative")
})
