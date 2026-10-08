test_that("per-group milestone survival and CI match survfit (all transforms)", {
  skip_if_not_installed("survival")

  set.seed(101)
  n_each <- 80
  time <- c(rexp(n_each, 0.12), rexp(n_each, 0.08))
  cens <- c(rexp(n_each, 0.04), rexp(n_each, 0.04))
  obs <- pmin(time, cens)
  status <- as.integer(time <= cens)
  group <- rep(c(0L, 1L), each = n_each)
  tau <- 8

  # conf.type used by survfit for each milestone_fast method's per-group CI.
  type_map <- c(wald = "plain", mover = "log", loglog = "log-log")

  for (m in names(type_map)) {
    res <- milestone_fast(obs, status, group, control = 0L, tau = tau, method = m)

    for (g in c(0L, 1L)) {
      idx <- group == g
      fit <- survival::survfit(
        survival::Surv(obs[idx], status[idx]) ~ 1,
        conf.type = type_map[[m]], conf.int = 0.95
      )
      s <- summary(fit, times = tau)
      key <- if (g == 0L) "control" else "treatment"

      expect_equal(unname(res$surv[key]), s$surv, tolerance = 1e-9,
                   info = paste(m, "surv group", g))
      expect_equal(unname(res$surv.lower[key]), s$lower, tolerance = 1e-6,
                   info = paste(m, "lower group", g))
      expect_equal(unname(res$surv.upper[key]), s$upper, tolerance = 1e-6,
                   info = paste(m, "upper group", g))
    }
  }
})

test_that("wald difference interval and statistic are computed correctly", {
  set.seed(202)
  n_each <- 60
  obs <- c(rexp(n_each, 0.15), rexp(n_each, 0.1))
  status <- rep(1L, 2 * n_each)
  group <- rep(c(0L, 1L), each = n_each)
  tau <- 6

  res <- milestone_fast(obs, status, group, control = 0L, tau = tau, method = "wald")

  z <- qnorm(0.975)
  v0 <- res$std.err["control"]^2
  v1 <- res$std.err["treatment"]^2
  se_diff <- sqrt(v0 + v1)
  d <- res$surv["treatment"] - res$surv["control"]

  expect_equal(res$diff, unname(d), tolerance = 1e-12)
  expect_equal(res$diff.lower, unname(d - z * se_diff), tolerance = 1e-12)
  expect_equal(res$diff.upper, unname(d + z * se_diff), tolerance = 1e-12)
  expect_equal(res$statistic, unname(d / se_diff), tolerance = 1e-12)
})

test_that("MOVER difference interval matches the one-sample recovery formula", {
  set.seed(303)
  n_each <- 70
  obs <- c(rexp(n_each, 0.12), rexp(n_each, 0.09))
  cens <- rexp(2 * n_each, 0.05)
  status <- as.integer(obs <= cens)
  obs <- pmin(obs, cens)
  group <- rep(c(0L, 1L), each = n_each)
  tau <- 7

  for (m in c("loglog", "mover")) {
    res <- milestone_fast(obs, status, group, control = 0L, tau = tau, method = m)

    s0 <- res$surv["control"]; s1 <- res$surv["treatment"]
    l0 <- res$surv.lower["control"]; u0 <- res$surv.upper["control"]
    l1 <- res$surv.lower["treatment"]; u1 <- res$surv.upper["treatment"]

    sigma_l <- (l1 - s1)^2 + (u0 - s0)^2
    sigma_u <- (u1 - s1)^2 + (l0 - s0)^2
    d <- s1 - s0

    expect_equal(res$diff.lower, unname(d - sqrt(sigma_l)), tolerance = 1e-12,
                 info = paste(m, "lower"))
    expect_equal(res$diff.upper, unname(d + sqrt(sigma_u)), tolerance = 1e-12,
                 info = paste(m, "upper"))
  }
})

test_that("loglog test statistic equals the variance-stabilised contrast", {
  set.seed(404)
  n_each <- 90
  obs <- c(rexp(n_each, 0.13), rexp(n_each, 0.08))
  status <- rep(1L, 2 * n_each)
  group <- rep(c(0L, 1L), each = n_each)
  tau <- 9

  res <- milestone_fast(obs, status, group, control = 0L, tau = tau, method = "loglog")

  s0 <- res$surv["control"]; s1 <- res$surv["treatment"]
  v0 <- res$std.err["control"]^2; v1 <- res$std.err["treatment"]^2
  g0 <- log(-log(s0)); g1 <- log(-log(s1))
  vg0 <- v0 / (s0 * log(s0))^2
  vg1 <- v1 / (s1 * log(s1))^2
  expected <- (g1 - g0) / sqrt(vg0 + vg1)

  expect_equal(res$statistic, unname(expected), tolerance = 1e-12)
})

test_that("presorted input gives identical results", {
  set.seed(505)
  n_each <- 50
  obs <- c(rexp(n_each, 0.1), rexp(n_each, 0.07))
  cens <- rexp(2 * n_each, 0.05)
  status <- as.integer(obs <= cens)
  obs <- pmin(obs, cens)
  group <- rep(c(0L, 1L), each = n_each)
  tau <- 8

  res_unsorted <- milestone_fast(obs, status, group, control = 0L, tau = tau, method = "loglog")

  ord <- order(obs)
  res_sorted <- milestone_fast(obs[ord], status[ord], group[ord], control = 0L,
                               tau = tau, method = "loglog", presorted = TRUE)

  expect_equal(res_unsorted$surv, res_sorted$surv, tolerance = 1e-12)
  expect_equal(res_unsorted$diff, res_sorted$diff, tolerance = 1e-12)
  expect_equal(res_unsorted$statistic, res_sorted$statistic, tolerance = 1e-12)
})

test_that("side argument controls the p-value direction", {
  set.seed(606)
  n_each <- 80
  obs <- c(rexp(n_each, 0.15), rexp(n_each, 0.08))
  status <- rep(1L, 2 * n_each)
  group <- rep(c(0L, 1L), each = n_each)
  tau <- 7

  two <- milestone_fast(obs, status, group, control = 0L, tau = tau, side = 2)
  one <- milestone_fast(obs, status, group, control = 0L, tau = tau, side = 1)

  # side = 1 reports the treatment-benefit one-sided p-value (upper tail for the
  # default Wald statistic); the two-sided p-value is twice the smaller tail.
  expect_equal(two$p.value, 2 * min(one$p.value, 1 - one$p.value),
               tolerance = 1e-10)
})

test_that("input validation catches malformed arguments", {
  obs <- c(1, 2, 3, 4)
  status <- c(1L, 0L, 1L, 1L)
  group <- c(0L, 0L, 1L, 1L)

  expect_error(milestone_fast(obs, status[-1], group, control = 0L, tau = 2),
               "same length")
  expect_error(milestone_fast(obs, c(1L, 2L, 1L, 1L), group, control = 0L, tau = 2),
               "0 .censored. and 1")
  expect_error(milestone_fast(obs, status, group, control = 0L, tau = -1),
               "positive")
  expect_error(milestone_fast(obs, status, rep(0L, 4), tau = 2),
               "two distinct")
  expect_error(milestone_fast(obs, status, group, control = 0L, tau = 2, conf.level = 1.2),
               "in .0, 1.")
  expect_error(milestone_fast(obs, status, group, control = 9, tau = 2),
               "must be one of")
})

test_that("print returns the object invisibly and produces output", {
  set.seed(707)
  n_each <- 40
  obs <- c(rexp(n_each, 0.1), rexp(n_each, 0.07))
  status <- rep(1L, 2 * n_each)
  group <- rep(c(0L, 1L), each = n_each)

  res <- milestone_fast(obs, status, group, control = 0L, tau = 8)
  expect_s3_class(res, "milestone_fast")
  expect_output(print(res), "Milestone survival")
  # capture.output keeps the printed report out of the test log
  utils::capture.output(vis <- withVisible(print(res)))
  expect_false(vis$visible)
  expect_identical(vis$value, res)
})

test_that("milestone_fast gives a zero standard error when a curve reaches zero", {
  # Control times 1, 2, 3 are all events, so S0(5) = 0; treatment has one
  # event at 4 with three at risk, so S1(5) = 2/3.
  # tau = 5 is beyond the control follow-up (largest time 3), so a warning
  # is expected as well.
  expect_warning(
    res <- milestone_fast(c(1, 2, 3, 4, 5, 6), c(1, 1, 1, 1, 0, 0),
                          c(0, 0, 0, 1, 1, 1), control = 0, tau = 5),
    "exceeds")
  expect_equal(unname(res$surv["control"]), 0)
  expect_equal(unname(res$std.err["control"]), 0)
  expect_equal(unname(res$surv["treatment"]), 2 / 3, tolerance = 1e-12)
})

test_that("milestone_fast rejects missing times and groups", {
  expect_error(milestone_fast(c(1, NA, 3, 4), c(1, 1, 0, 1), c(0, 0, 1, 1),
                              control = 0, tau = 2), "missing")
  expect_error(milestone_fast(c(1, 2, 3, 4), c(1, 1, 0, 1), c(0, NA, 1, 1),
                              control = 0, tau = 2), "missing")
})

test_that("milestone_fast: an estimate of 1 gives a degenerate one-sample interval", {
  set.seed(3)
  # No control event before tau = 5, so the control estimate is 1.
  time  <- c(runif(40, 10, 20), rexp(40, 0.1))
  event <- rep(1L, 80)
  group <- rep(0:1, each = 40)
  d <- data.frame(sim = 1L, group = group + 1L, accrual_time = 0,
                  tte = time, event = event)
  for (m in c("mover", "loglog")) {
    r <- milestone_fast(time, event, group, control = 0, tau = 5, method = m)
    expect_equal(unname(r$surv[1]), 1)
    expect_equal(unname(c(r$surv.lower[1], r$surv.upper[1])), c(1, 1))
    expect_true(is.finite(r$diff.lower) && is.finite(r$diff.upper))
    expect_true(r$diff.lower <= r$diff && r$diff <= r$diff.upper)
    a <- analysis_fast(d, control = 1, time.looks = 1000, stat = "milestone",
                       tau = 5, ms.method = m)
    expect_equal(a$milestone.diff.lower, r$diff.lower)
    expect_equal(a$milestone.diff.upper, r$diff.upper)
  }
  expect_true(is.finite(milestone_fast(time, event, group, control = 0,
                                       tau = 5, method = "mover")$statistic))
  expect_true(is.na(milestone_fast(time, event, group, control = 0,
                                   tau = 5, method = "loglog")$statistic))
})

test_that("milestone_fast: event values are validated before conversion", {
  expect_error(milestone_fast(1:4, c(1, 0.7, 0, 1), c(0, 0, 1, 1),
                              control = 0, tau = 2), "only 0")
  expect_error(milestone_fast(c(-1, 2, 3, 4), c(1, 1, 0, 1), c(0, 0, 1, 1),
                              control = 0, tau = 2), "non-negative")
})

test_that("milestone_fast: presorted = TRUE checks the order", {
  expect_error(milestone_fast(c(3, 1, 2, 4), c(1, 1, 0, 1), c(0, 1, 0, 1),
                              control = 0, tau = 2, presorted = TRUE),
               "presorted = FALSE")
})
