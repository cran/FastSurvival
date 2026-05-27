test_that("coxph_fast returns named vector of length 5", {
  set.seed(1)
  n     <- 100
  time  <- c(rexp(n / 2, 0.1), rexp(n / 2, 0.15))
  event <- rep(1L, n)
  group <- rep(c(1L, 2L), each = n / 2)

  res <- coxph_fast(time, event, group, control = 1)
  expect_length(res, 5L)
  expect_named(res, c("coef", "exp(coef)", "se(coef)", "lower .95", "upper .95"))
})

test_that("coxph_fast exp(coef) is positive", {
  set.seed(2)
  n     <- 100
  time  <- c(rexp(n / 2, 0.1), rexp(n / 2, 0.15))
  event <- rep(1L, n)
  group <- rep(c(1L, 2L), each = n / 2)

  res <- coxph_fast(time, event, group, control = 1)
  expect_gt(unname(res["exp(coef)"]), 0)
})

test_that("coxph_fast CI lower <= exp(coef) <= upper", {
  set.seed(3)
  n     <- 100
  time  <- c(rexp(n / 2, 0.1), rexp(n / 2, 0.15))
  event <- rep(1L, n)
  group <- rep(c(1L, 2L), each = n / 2)

  res <- coxph_fast(time, event, group, control = 1)
  expect_lte(unname(res["lower .95"]), unname(res["exp(coef)"]) + 1e-10)
  expect_gte(unname(res["upper .95"]), unname(res["exp(coef)"]) - 1e-10)
})

test_that("coxph_fast coef = log(exp(coef))", {
  set.seed(4)
  n     <- 100
  time  <- c(rexp(n / 2, 0.1), rexp(n / 2, 0.15))
  event <- rep(1L, n)
  group <- rep(c(1L, 2L), each = n / 2)

  res <- coxph_fast(time, event, group, control = 1)
  expect_equal(unname(res["coef"]), log(unname(res["exp(coef)"])),
               tolerance = 1e-12)
})

test_that("coxph_fast presorted=TRUE and presorted=FALSE agree", {
  set.seed(5)
  n     <- 100
  time  <- c(rexp(n / 2, 0.1), rexp(n / 2, 0.15))
  event <- rep(1L, n)
  group <- rep(c(1L, 2L), each = n / 2)
  ord   <- order(time)

  res_pre <- coxph_fast(time[ord], event[ord], group[ord],
                        control = 1, presorted = TRUE)
  res_uns <- coxph_fast(time, event, group,
                        control = 1, presorted = FALSE)
  expect_equal(unname(res_pre), unname(res_uns), tolerance = 1e-12)
})

test_that("coxph_fast conf.int argument changes CI label names", {
  set.seed(6)
  n     <- 100
  time  <- c(rexp(n / 2, 0.1), rexp(n / 2, 0.15))
  event <- rep(1L, n)
  group <- rep(c(1L, 2L), each = n / 2)

  res90 <- coxph_fast(time, event, group, control = 1, conf.int = 0.90)
  expect_named(res90, c("coef", "exp(coef)", "se(coef)", "lower .90", "upper .90"))
})

test_that("coxph_fast handles factor group correctly", {
  set.seed(7)
  n     <- 100
  time  <- c(rexp(n / 2, 0.1), rexp(n / 2, 0.15))
  event <- rep(1L, n)
  group <- factor(rep(c("control", "treatment"), each = n / 2))

  res <- coxph_fast(time, event, group, control = "control")
  expect_gt(unname(res["exp(coef)"]), 0)
  expect_false(any(is.na(res)))
})

test_that("coxph_fast agrees with coxph (no ties, tolerance 1e-4)", {
  skip_if_not_installed("survival")
  set.seed(8)
  n     <- 200
  time  <- c(rexp(n / 2, 0.1), rexp(n / 2, 0.15))
  event <- rep(1L, n)
  group <- rep(c(1L, 2L), each = n / 2)

  res <- coxph_fast(time, event, group, control = 1)

  fit <- survival::coxph(
    survival::Surv(time, event) ~ I(group != 1),
    ties    = "breslow",
    control = survival::coxph.control(
      eps        = 1e-12,
      toler.chol = .Machine$double.eps^0.875
    )
  )
  hr_cox <- unname(exp(stats::coef(fit)[1L]))

  expect_equal(unname(res["exp(coef)"]), hr_cox, tolerance = 1e-4)
})

test_that("coxph_fast returns all-NA when no events", {
  time  <- c(1, 2, 3, 4)
  event <- c(0, 0, 0, 0)
  group <- c(1, 1, 2, 2)

  res <- coxph_fast(time, event, group, control = 1)
  expect_true(all(is.na(res)))
})

test_that("coxph_fast stops when input lengths differ", {
  expect_error(
    coxph_fast(c(1, 2, 3), c(1, 0), c(1, 1, 2), control = 1),
    "same length"
  )
})
