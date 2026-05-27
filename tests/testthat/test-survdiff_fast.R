test_that("survdiff_fast returns a single numeric value", {
  set.seed(1)
  n     <- 100
  time  <- c(rexp(n / 2, 0.1), rexp(n / 2, 0.15))
  event <- rep(1L, n)
  group <- rep(c(1L, 2L), each = n / 2)

  res <- survdiff_fast(time, event, group, control = 1, side = 2)
  expect_length(res, 1L)
  expect_true(is.numeric(res))
})

test_that("survdiff_fast chi-square is non-negative", {
  set.seed(2)
  n     <- 100
  time  <- c(rexp(n / 2, 0.1), rexp(n / 2, 0.15))
  event <- rep(1L, n)
  group <- rep(c(1L, 2L), each = n / 2)

  res <- survdiff_fast(time, event, group, control = 1, side = 2)
  expect_gte(as.numeric(res), 0)
})

test_that("survdiff_fast side=1 and side=2 are consistent (Z^2 = chi-sq)", {
  set.seed(3)
  n     <- 100
  time  <- c(rexp(n / 2, 0.1), rexp(n / 2, 0.15))
  event <- rep(1L, n)
  group <- rep(c(1L, 2L), each = n / 2)

  z   <- survdiff_fast(time, event, group, control = 1, side = 1)
  chi <- survdiff_fast(time, event, group, control = 1, side = 2)
  expect_equal(as.numeric(z) ^ 2, as.numeric(chi), tolerance = 1e-10)
})

test_that("survdiff_fast presorted=TRUE and presorted=FALSE agree", {
  set.seed(4)
  n     <- 100
  time  <- c(rexp(n / 2, 0.1), rexp(n / 2, 0.15))
  event <- rep(1L, n)
  group <- rep(c(1L, 2L), each = n / 2)
  ord   <- order(time)

  res_pre <- survdiff_fast(time[ord], event[ord], group[ord],
                           control = 1, side = 2, presorted = TRUE)
  res_uns <- survdiff_fast(time, event, group,
                           control = 1, side = 2, presorted = FALSE)
  expect_equal(as.numeric(res_pre), as.numeric(res_uns), tolerance = 1e-10)
})

test_that("survdiff_fast result is symmetric under group relabeling", {
  # Swapping control and treatment negates the Z-score but preserves chi-square
  set.seed(9)
  n     <- 100
  time  <- c(rexp(n / 2, 0.1), rexp(n / 2, 0.15))
  event <- rep(1L, n)
  group <- rep(c(1L, 2L), each = n / 2)

  z_12 <- survdiff_fast(time, event, group, control = 1, side = 1)
  z_21 <- survdiff_fast(time, event, group, control = 2, side = 1)
  chi  <- survdiff_fast(time, event, group, control = 1, side = 2)

  expect_equal(as.numeric(z_12), -as.numeric(z_21), tolerance = 1e-10)
  expect_equal(as.numeric(z_12) ^ 2, as.numeric(chi), tolerance = 1e-10)
})

test_that("survdiff_fast agrees with survival::survdiff chi-square", {
  skip_if_not_installed("survival")
  set.seed(5)
  n     <- 200
  time  <- c(rexp(n / 2, 0.1), rexp(n / 2, 0.15))
  event <- rbinom(n, 1, 0.8)
  group <- rep(c(1L, 2L), each = n / 2)

  res <- survdiff_fast(time, event, group, control = 1, side = 2)
  ref <- survival::survdiff(survival::Surv(time, event) ~ group)

  expect_equal(as.numeric(res), ref$chisq, tolerance = 1e-6)
})

test_that("survdiff_fast stops when no events observed", {
  time  <- c(1, 2, 3, 4)
  event <- c(0, 0, 0, 0)
  group <- c(1, 1, 2, 2)

  expect_error(
    survdiff_fast(time, event, group, control = 1, side = 2),
    "No events"
  )
})

test_that("survdiff_fast stops when input lengths differ", {
  expect_error(
    survdiff_fast(c(1, 2, 3), c(1, 0), c(1, 1, 2), control = 1, side = 2),
    "same length"
  )
})
