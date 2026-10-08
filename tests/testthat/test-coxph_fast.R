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

test_that("coxph_fast conf.level argument changes CI label names", {
  set.seed(6)
  n     <- 100
  time  <- c(rexp(n / 2, 0.1), rexp(n / 2, 0.15))
  event <- rep(1L, n)
  group <- rep(c(1L, 2L), each = n / 2)

  res90 <- coxph_fast(time, event, group, control = 1, conf.level = 0.90)
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

test_that("coxph_fast validates event coding, group levels, and conf.level", {
  tt <- c(1, 2, 3, 4, 5, 6)
  gg <- c(1, 1, 1, 2, 2, 2)
  expect_error(coxph_fast(tt, c(1, 2, 1, 2, 1, 2), gg, control = 1),
               "coded as 0")
  expect_error(coxph_fast(tt, rep(1, 6), c(1, 1, 2, 2, 3, 3), control = 1),
               "two distinct")
  expect_error(coxph_fast(tt, rep(1, 6), gg, control = 3), "control")
  expect_error(coxph_fast(tt, rep(1, 6), gg, control = 1, conf.level = 95),
               "conf.level")
})

# ------------------------------------------------------------------ #
#  Stratified estimator
# ------------------------------------------------------------------ #

test_that("stratified core reproduces a hand calculation", {
  # Stratum 1: (t, group, event) = (1, C, 1), (2, T, 1), (3, C, 0), (4, T, 1)
  # Stratum 2: (1, T, 1), (2, C, 1), (3, T, 1)
  # Risk sets restart in each stratum, so O_T = 4, O_C = 2,
  # E_T = 2/4 + 2/3 + 1 + 2/3 + 1/2 + 1 = 13/3 and
  # E_C = 2/4 + 1/3 + 0 + 1/3 + 1/2 + 0 = 5/3, giving the Pike anchor
  # theta_0 = (4 * 5/3) / (2 * 13/3) = 10/13. U_0 = -62/759 and
  # I_0 = 558220/576081 follow from p_k = n_T theta_0 / (n_C + n_T theta_0).
  res <- pihe_core_strat(c(1, 2, 3, 4, 1, 2, 3), c(1L, 1L, 0L, 1L, 1L, 1L, 1L),
                         c(0L, 1L, 0L, 1L, 1L, 0L, 1L),
                         c(1L, 1L, 1L, 1L, 2L, 2L, 2L))
  expect_equal(res[1], 10 / 13, tolerance = 1e-12)
  expect_equal(res[2], -62 / 759, tolerance = 1e-12)
  expect_equal(res[3], 558220 / 576081, tolerance = 1e-12)
})

test_that("a single stratum gives the unstratified estimate", {
  set.seed(401)
  n  <- 200
  tt <- rexp(n, 0.1)
  ee <- rbinom(n, 1, 0.8)
  gg <- rep(1:2, each = n / 2)
  un <- coxph_fast(tt, ee, gg, control = 1)
  st <- coxph_fast(tt, ee, gg, control = 1, strata = rep("all", n))
  expect_equal(as.numeric(st), as.numeric(un), tolerance = 1e-12)
  expect_equal(attr(st, "strata"), 1L)
})

test_that("stratified coxph_fast agrees with survival::coxph(strata())", {
  skip_if_not_installed("survival")
  set.seed(402)
  n_s  <- 150
  base <- c(0.05, 0.10, 0.20, 0.40)
  s    <- rep(seq_along(base), each = n_s)
  g    <- rep(0:1, times = length(base) * n_s / 2)
  lam  <- base[s] * ifelse(g == 1, 0.7, 1)
  t0   <- rexp(length(s), lam)
  cc   <- rexp(length(s), 0.05)
  tt   <- round(pmin(t0, cc), 1)       # rounding creates tied times
  ee   <- as.integer(t0 <= cc)
  fast <- coxph_fast(tt, ee, g, control = 0, strata = s)
  ref  <- survival::coxph(survival::Surv(tt, ee) ~ g + survival::strata(s),
                          ties = "breslow")
  expect_lt(abs(unname(fast["coef"]) - unname(stats::coef(ref))), 1e-3)
  expect_equal(unname(fast["se(coef)"]), unname(sqrt(ref$var[1, 1])),
               tolerance = 1e-2)
  expect_equal(attr(fast, "strata"), 4L)
})

test_that("stratified coxph_fast does not depend on the row order", {
  set.seed(403)
  n  <- 240
  s  <- sample(c("a", "b", "c"), n, replace = TRUE)
  tt <- rexp(n, ifelse(s == "a", 0.05, 0.2))
  ee <- rbinom(n, 1, 0.8)
  gg <- sample(1:2, n, replace = TRUE)
  ref <- coxph_fast(tt, ee, gg, control = 1, strata = s)
  shuf <- sample(n)
  expect_equal(as.numeric(coxph_fast(tt[shuf], ee[shuf], gg[shuf],
                                     control = 1, strata = s[shuf])),
               as.numeric(ref), tolerance = 1e-12)
  ord <- order(s, tt)
  expect_equal(as.numeric(coxph_fast(tt[ord], ee[ord], gg[ord], control = 1,
                                     strata = s[ord], presorted = TRUE)),
               as.numeric(ref), tolerance = 1e-12)
  expect_output(print(ref), "Stratified Pike-Halley estimator")
})

test_that("stratified coxph_fast validates strata", {
  tt <- c(1, 2, 3, 4, 5, 6)
  gg <- c(1, 1, 1, 2, 2, 2)
  ee <- rep(1, 6)
  expect_error(coxph_fast(tt, ee, gg, control = 1, strata = c(1, 2)),
               "same length")
  expect_error(coxph_fast(tt, ee, gg, control = 1,
                          strata = c(1, NA, 1, 2, 2, 2)), "missing")
  expect_error(coxph_fast(tt, ee, gg, control = 1, presorted = TRUE,
                          strata = c(1, 2, 1, 2, 1, 2)), "contiguous")
})

test_that("coxph_fast: presorted = TRUE checks the order", {
  time  <- c(3, 1, 2, 4, 6, 5)
  event <- c(1, 1, 0, 1, 1, 0)
  group <- c(0, 1, 0, 1, 0, 1)
  expect_error(coxph_fast(time, event, group, control = 0, presorted = TRUE), "presorted = FALSE")
  expect_error(coxph_fast(time, event, group, control = 0, presorted = TRUE, strata = c(1, 1, 1, 2, 2, 2)), "presorted = FALSE")
  expect_error(coxph_fast(c(1, 3, 2, 4, 5, 6), event, group, control = 0, presorted = TRUE, strata = c(1, 2, 1, 2, 1, 2)), "presorted = FALSE")
})
