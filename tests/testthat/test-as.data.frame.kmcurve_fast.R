test_that("as.data.frame.kmcurve_fast matches a hand calculation", {
  # Hand calculation with exact fractions: S = 5/6, 2/3, 4/9, 4/9, 0 and
  # Greenwood sums 1/30, 1/12, 1/4, 1/4, 1/4
  # (no term when all subjects at risk fail)
  fit <- kmcurve_fast(c(1, 2, 2, 3, 4, 5), c(1, 1, 0, 1, 0, 1))
  tab <- as.data.frame(fit)
  expect_identical(names(tab), c("group", "time", "n.risk", "n.event",
                                 "n.censor", "surv", "std.err", "lower",
                                 "upper"))
  expect_identical(tab$group, rep("All", 6L))
  expect_equal(tab$time, c(0, 1, 2, 3, 4, 5))
  expect_equal(tab$n.risk, c(6, 6, 5, 3, 2, 1))
  expect_equal(tab$n.event, c(0, 1, 1, 1, 0, 1))
  expect_equal(tab$n.censor, c(0, 0, 1, 0, 1, 0))
  expect_equal(tab$surv, c(1, 5 / 6, 2 / 3, 4 / 9, 4 / 9, 0))
  expect_equal(tab$std.err,
               c(0, 0.152145154862546, 0.192450089729875, 0.222222222222222,
                 0.222222222222222, 0), tolerance = 1e-12)
  expect_equal(tab$lower[2:5],
               c(0.582654795477, 0.378606460887, 0.166807936628,
                 0.166807936628), tolerance = 1e-10)
  expect_equal(tab$upper, c(1, 1, 1, 1, 1, 0))
  expect_identical(rownames(as.data.frame(fit, row.names = letters[1:6])),
                   letters[1:6])
})

test_that("as.data.frame.kmcurve_fast agrees with the curves of the object", {
  set.seed(95)
  n <- 60
  time  <- round(rexp(2 * n, 0.1), 1)
  event <- rbinom(2 * n, 1, 0.7)
  group <- rep(c("A", "B"), each = n)
  fit <- kmcurve_fast(time, event, group, control = "B")
  tab <- as.data.frame(fit, conf.level = 0.9, conf.type = "plain")
  # The control group comes first
  expect_identical(unique(tab$group), c("B", "A"))
  for (role in c("control", "treat")) {
    km <- fit$km[[role]]
    tg <- tab[tab$group == fit$labels[[role]], ]
    ev <- tg$n.event > 0
    expect_equal(tg$time[ev], km$te)
    expect_equal(tg$surv[ev], km$surv)
    expect_equal(tg$std.err[ev], km$se)
    ci <- km_ci(km$surv, km$se, 0.9, "plain")
    expect_equal(tg$lower[ev], ci$lower)
    expect_equal(tg$upper[ev], ci$upper)
    expect_equal(tg$n.risk[1L], km$n)
    expect_equal(sum(tg$n.event), km$n_event)
  }
  expect_error(as.data.frame(fit, conf.level = 1), "conf.level")
  expect_error(as.data.frame(fit, conf.type = "arcsin"), "should be one of")
})

test_that("as.data.frame.kmcurve_fast matches survival::survfit", {
  set.seed(96)
  n  <- 100
  t0 <- rexp(n, log(2) / 12)
  cc <- runif(n, 0, 30)
  time  <- c(round(pmin(t0, cc), 1), 5, 5, 5)
  event <- c(as.integer(t0 <= cc), 1, 0, 1)
  tab <- as.data.frame(kmcurve_fast(time, event), conf.type = "log-log")
  # Report whether the external comparison below is run or skipped.
  message("survival available for the as.data.frame.kmcurve_fast ",
          "comparison: ", requireNamespace("survival", quietly = TRUE))
  skip_if_not_installed("survival")
  ref <- survival::survfit(survival::Surv(time, event) ~ 1,
                           conf.type = "log-log")
  sm  <- summary(ref, censored = TRUE)
  tb  <- tab[-1L, ]
  expect_equal(tb$time, sm$time)
  expect_equal(tb$n.risk, sm$n.risk)
  expect_equal(tb$n.event, sm$n.event)
  expect_equal(tb$n.censor, sm$n.censor)
  expect_equal(tb$surv, sm$surv, tolerance = 1e-10)
  inner <- tb$surv > 0 & tb$surv < 1
  expect_equal(tb$std.err[inner], sm$std.err[inner], tolerance = 1e-8)
  expect_equal(tb$lower[inner], sm$lower[inner], tolerance = 1e-6)
  expect_equal(tb$upper[inner], sm$upper[inner], tolerance = 1e-6)
})
