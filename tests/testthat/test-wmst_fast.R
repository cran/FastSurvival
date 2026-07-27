# Reference window mean survival time and its variance from survfit output,
# using the same Greenwood form as wmst_core. Provides an independent
# validation axis. tt, ss: time and event (1 = event, 0 = censor).
wmst_reference <- function(tt, ss, tau1, tau2) {
  ft <- survival::survfit(survival::Surv(tt, ss) ~ 1)
  keep <- ft$n.event > 0 & ft$time <= tau2
  e <- ft$time[keep]
  s_after <- ft$surv[keep]
  y <- ft$n.risk[keep]
  d <- ft$n.event[keep]
  m <- length(e)

  left_t <- c(0, e)
  right_t <- c(e, tau2)
  val <- c(1, s_after)
  lo <- pmax(left_t, tau1)
  hi <- pmin(right_t, tau2)
  wmst <- sum(val * pmax(0, hi - lo))

  width_to_tau2 <- pmax(0, pmin(right_t, tau2) - left_t)
  area_i <- val * width_to_tau2
  tail_cum <- rev(cumsum(rev(area_i)))
  t_k <- if (m >= 1) tail_cum[2:(m + 1)] else numeric(0)
  rw <- ifelse(e >= tau1, t_k, wmst)
  vv <- ifelse(y - d == 0, 0, d / (y * (y - d)))
  variance <- sum(rw^2 * vv)

  list(wmst = wmst, variance = variance)
}

test_that("wmst_fast returns the expected structure for two groups", {
  set.seed(1)
  n <- 200
  g <- rep(0:1, each = n / 2)
  tt <- c(rexp(n / 2, log(2) / 12), rexp(n / 2, log(2) / 16))
  cc <- rexp(n, rate = 0.02)
  time <- pmin(tt, cc)
  event <- as.integer(tt <= cc)

  res <- wmst_fast(time, event, group = g, control = 0, tau1 = 2, tau2 = 18)
  expect_s3_class(res, "wmst_fast")
  expect_true(all(c("wmst.control", "wmst.treatment", "diff",
                    "se.diff", "z", "chisq", "p") %in% names(res)))
  expect_equal(unname(res["chisq"]), unname(res["z"])^2, tolerance = 1e-8)
  expect_true(is.finite(res["p"]) && res["p"] >= 0 && res["p"] <= 1)
})

test_that("wmst_fast with tau1 = 0 matches survRM2::rmst2", {
  skip_if_not_installed("survRM2")
  out <- tryCatch({
    set.seed(21)
    n_per <- 250
    g <- rep(0:1, each = n_per)
    tt <- c(rexp(n_per, log(2) / 12), rexp(n_per, log(2) / 16))
    cc <- rexp(2 * n_per, rate = 0.02)
    time <- pmin(tt, cc)
    event <- as.integer(tt <= cc)
    tau <- min(max(time[g == 0]), max(time[g == 1]))
    rm2 <- survRM2::rmst2(time = time, status = event, arm = g, tau = tau)
    diff_ref <- as.numeric(rm2$RMST.arm1$rmst["Est."] - rm2$RMST.arm0$rmst["Est."])
    se_ref <- sqrt(rm2$RMST.arm1$rmst.var + rm2$RMST.arm0$rmst.var)
    fast <- wmst_fast(time, event, group = g, control = 0, tau1 = 0, tau2 = tau)
    list(diff_ref = diff_ref, se_ref = se_ref, fast = fast)
  }, error = function(e) NULL)
  skip_if(is.null(out), "survRM2 comparison unavailable")
  expect_equal(unname(out$fast["diff"]), out$diff_ref, tolerance = 1e-6)
  expect_equal(unname(out$fast["se.diff"]), out$se_ref, tolerance = 1e-6)
})

test_that("wmst_fast matches the survfit-based reference over a window", {
  skip_if_not_installed("survival")
  out <- tryCatch({
    set.seed(31)
    n_per <- 250
    g <- rep(0:1, each = n_per)
    tt <- c(rexp(n_per, log(2) / 12), rexp(n_per, log(2) / 16))
    cc <- rexp(2 * n_per, rate = 0.02)
    time <- pmin(tt, cc)
    event <- as.integer(tt <= cc)
    tau1 <- 3
    tau2 <- min(max(time[g == 0]), max(time[g == 1]))
    r0 <- wmst_reference(time[g == 0], event[g == 0], tau1, tau2)
    r1 <- wmst_reference(time[g == 1], event[g == 1], tau1, tau2)
    fast <- wmst_fast(time, event, group = g, control = 0,
                      tau1 = tau1, tau2 = tau2)
    list(r0 = r0, r1 = r1, fast = fast)
  }, error = function(e) NULL)
  skip_if(is.null(out), "survival reference unavailable")
  expect_equal(unname(out$fast["wmst.control"]), out$r0$wmst, tolerance = 1e-8)
  expect_equal(unname(out$fast["wmst.treatment"]), out$r1$wmst, tolerance = 1e-8)
  expect_equal(unname(out$fast["se.diff"]),
               sqrt(out$r0$variance + out$r1$variance), tolerance = 1e-8)
})

test_that("a treatment benefit gives a positive difference and z", {
  set.seed(7)
  n_per <- 300
  g <- rep(0:1, each = n_per)
  tt <- c(rexp(n_per, log(2) / 10), rexp(n_per, log(2) / 18))
  cc <- rexp(2 * n_per, rate = 0.02)
  time <- pmin(tt, cc)
  event <- as.integer(tt <= cc)

  res <- wmst_fast(time, event, group = g, control = 0, tau1 = 0, tau2 = 14,
                   side = 1)
  expect_gt(unname(res["diff"]), 0)
  expect_gt(unname(res["z"]), 0)
})

test_that("one-sided and two-sided p-values are consistent with z", {
  set.seed(8)
  n_per <- 250
  g <- rep(0:1, each = n_per)
  tt <- c(rexp(n_per, log(2) / 10), rexp(n_per, log(2) / 16))
  cc <- rexp(2 * n_per, rate = 0.02)
  time <- pmin(tt, cc)
  event <- as.integer(tt <= cc)

  res2 <- wmst_fast(time, event, group = g, control = 0, tau2 = 14, side = 2)
  res1 <- wmst_fast(time, event, group = g, control = 0, tau2 = 14, side = 1)
  z <- unname(res2["z"])
  expect_equal(unname(res1["p"]), pnorm(z, lower.tail = FALSE), tolerance = 1e-8)
  expect_equal(unname(res2["p"]), 2 * pnorm(-abs(z)), tolerance = 1e-8)
})

test_that("single-group mode returns a WMST and confidence interval", {
  set.seed(5)
  n <- 200
  tt <- rexp(n, 0.1)
  cc <- rexp(n, 0.02)
  time <- pmin(tt, cc)
  event <- as.integer(tt <= cc)

  res <- wmst_fast(time, event, tau1 = 2, tau2 = 12)
  expect_s3_class(res, "wmst_fast")
  expect_true(all(c("wmst", "se", "lower", "upper") %in% names(res)))
  expect_false("p" %in% names(res))
  expect_true(unname(res["lower"]) <= unname(res["wmst"]))
  expect_true(unname(res["upper"]) >= unname(res["wmst"]))
})

test_that("input validation works", {
  expect_error(wmst_fast(1:5, c(0, 1, 0, 1)), "same length")
  expect_error(wmst_fast(1:5, rep(2L, 5)), "0")
  expect_error(
    wmst_fast(1:6, rep(0:1, 3), group = rep(1:3, 2), control = 1), "two distinct"
  )
  expect_error(
    wmst_fast(1:6, rep(0:1, 3), group = rep(0:1, 3)), "control must be specified"
  )
  expect_error(
    wmst_fast(1:6, rep(0:1, 3), group = rep(0:1, 3), control = 0,
              tau1 = 5, tau2 = 3), "tau2 must be greater"
  )
  expect_error(
    wmst_fast(1:6, rep(0:1, 3), group = rep(0:1, 3), control = 0, tau1 = -1),
    "tau1 must be non-negative"
  )
})

test_that("type I error is approximately controlled under the null", {
  set.seed(42)
  nsim <- 200L
  reject <- 0L
  for (s in seq_len(nsim)) {
    n_per <- 150
    g <- rep(0:1, each = n_per)
    tt <- rexp(2 * n_per, rate = log(2) / 12)
    cc <- rexp(2 * n_per, rate = 0.02)
    time <- pmin(tt, cc)
    event <- as.integer(tt <= cc)
    res <- wmst_fast(time, event, group = g, control = 0, tau2 = 14, side = 2)
    if (is.finite(res["p"]) && res["p"] < 0.05) reject <- reject + 1L
  }
  expect_lt(reject / nsim, 0.12)
})
