test_that("rmst_fast single-group returns named vector of length 4", {
  set.seed(1)
  t_raw <- rexp(100, rate = 1 / 10)
  e_raw <- rbinom(100, 1, 0.7)

  res <- rmst_fast(t_raw, e_raw, tau = 10)
  expect_length(res, 4L)
  expect_named(res, c("rmst", "std.err", "lower", "upper"))
})

test_that("rmst_fast single-group RMST is in [0, tau]", {
  set.seed(2)
  t_raw <- rexp(100, rate = 1 / 10)
  e_raw <- rbinom(100, 1, 0.7)

  res <- rmst_fast(t_raw, e_raw, tau = 10)
  expect_gte(unname(res["rmst"]), 0)
  expect_lte(unname(res["rmst"]), 10)
})

test_that("rmst_fast two-group returns difference and ratio contrasts", {
  set.seed(3)
  n     <- 200
  time  <- c(rexp(n, 0.10), rexp(n, 0.07))
  event <- rbinom(2 * n, 1, 0.8)
  group <- rep(0:1, each = n)

  res <- rmst_fast(time, event, group = group, control = 0, tau = 10)
  expect_true(all(c("rmst.ctrl", "rmst.trt", "diff", "se.diff",
                    "diff.lower", "diff.upper", "z.diff", "p.diff",
                    "ratio", "ratio.lower", "ratio.upper",
                    "z.ratio", "p.ratio") %in% names(res)))
})

test_that("rmst_fast diff equals treatment minus control", {
  set.seed(4)
  n     <- 200
  time  <- c(rexp(n, 0.10), rexp(n, 0.07))
  event <- rbinom(2 * n, 1, 0.8)
  group <- rep(0:1, each = n)

  res <- rmst_fast(time, event, group = group, control = 0, tau = 10)
  expect_equal(unname(res["diff"]),
               unname(res["rmst.trt"] - res["rmst.ctrl"]),
               tolerance = 1e-10)
})

test_that("rmst_fast two-sided p.diff equals 2 * pnorm(-|z.diff|)", {
  set.seed(5)
  n     <- 200
  time  <- c(rexp(n, 0.10), rexp(n, 0.07))
  event <- rbinom(2 * n, 1, 0.8)
  group <- rep(0:1, each = n)

  res <- rmst_fast(time, event, group = group, control = 0, tau = 10)
  expect_equal(unname(res["p.diff"]),
               2 * pnorm(-abs(unname(res["z.diff"]))),
               tolerance = 1e-12)
})

test_that("rmst_fast presorted=TRUE and presorted=FALSE agree", {
  set.seed(6)
  n     <- 200
  time  <- c(rexp(n, 0.10), rexp(n, 0.07))
  event <- rbinom(2 * n, 1, 0.8)
  group <- rep(0:1, each = n)
  ord   <- order(time)

  res_uns <- rmst_fast(time, event, group = group, control = 0, tau = 10)
  res_pre <- rmst_fast(time[ord], event[ord], group = group[ord],
                       control = 0, tau = 10, presorted = TRUE)
  expect_equal(as.numeric(res_uns), as.numeric(res_pre), tolerance = 1e-10)
})

test_that("rmst_fast single-group agrees with survival::survfit rmean", {
  skip_if_not_installed("survival")
  set.seed(7)
  t_raw <- rexp(150, rate = 1 / 10)
  e_raw <- rbinom(150, 1, 0.75)
  tau   <- 10

  res <- rmst_fast(t_raw, e_raw, tau = tau)

  fit <- survival::survfit(survival::Surv(t_raw, e_raw) ~ 1)
  tab <- summary(fit, rmean = tau)$table
  rmean <- unname(tab["rmean"])

  expect_equal(unname(res["rmst"]), rmean, tolerance = 1e-6)
})

test_that("rmst_fast two-group contrasts agree with survRM2", {
  skip_if_not_installed("survRM2")
  set.seed(8)
  n     <- 200
  time  <- c(rexp(n, 0.10), rexp(n, 0.07))
  event <- rbinom(2 * n, 1, 0.8)
  arm   <- rep(0:1, each = n)
  tau   <- 10

  res <- rmst_fast(time, event, group = arm, control = 0, tau = tau)
  ref <- survRM2::rmst2(time = time, status = event, arm = arm, tau = tau)
  ur  <- ref$unadjusted.result
  # ur rows: 1 = difference (arm1 - arm0), 2 = ratio (arm1 / arm0)
  expect_equal(unname(res["diff"]),       ur[1L, 1L], tolerance = 1e-6)
  expect_equal(unname(res["diff.lower"]), ur[1L, 2L], tolerance = 1e-6)
  expect_equal(unname(res["diff.upper"]), ur[1L, 3L], tolerance = 1e-6)
  expect_equal(unname(res["p.diff"]),     ur[1L, 4L], tolerance = 1e-6)
  expect_equal(unname(res["ratio"]),      ur[2L, 1L], tolerance = 1e-6)
})

test_that("rmst_fast single-group returns NA with class when n = 0", {
  res <- rmst_fast(numeric(0), integer(0), tau = 10)
  expect_s3_class(res, "rmst_fast")
  expect_true(all(is.na(res)))
})

test_that("rmst_fast validates event coding and the control label", {
  tt <- c(1, 2, 3, 4, 5, 6)
  gg <- c(0, 0, 0, 1, 1, 1)
  expect_error(rmst_fast(tt, c(1, 2, 1, 2, 1, 2), tau = 4), "coded as 0")
  expect_error(rmst_fast(tt, rep(1, 6), group = gg, control = 7, tau = 4),
               "control")
})

test_that("rmst_fast warns when tau exceeds the observed follow-up", {
  set.seed(71)
  tt <- c(rexp(50, 0.10), rexp(50, 0.07))
  ee <- rep(1L, 100)
  gg <- rep(0:1, each = 50)
  lim <- min(max(tt[gg == 0]), max(tt[gg == 1]))
  expect_warning(rmst_fast(tt, ee, group = gg, control = 0, tau = lim + 1),
                 "exceeds")
  expect_warning(rmst_fast(tt, ee, group = gg, control = 0, tau = lim), NA)
  expect_warning(rmst_fast(tt, ee, tau = max(tt) + 1), "exceeds")
})

test_that("rmst_fast: presorted = TRUE checks the order", {
  time  <- c(3, 1, 2, 4, 6, 5)
  event <- c(1, 1, 0, 1, 1, 0)
  group <- c(0, 1, 0, 1, 0, 1)
  expect_error(rmst_fast(time, event, tau = 3, presorted = TRUE), "presorted = FALSE")
  expect_error(rmst_fast(time, event, group, control = 0, tau = 3, presorted = TRUE), "presorted = FALSE")
})
