# Reference implementation of the Pepe-Fleming weighted Kaplan-Meier statistic,
# built from survfit and approxfun. This mirrors the algorithm of
# nphsim::wkm.Stat and provides an independent validation axis for wkm_core.
# event: 1 = event, 0 = censoring; grp: 1 = treatment, 0 = control.
# The selected weight enters the numerator and the weighted area; the variance
# always divides by the Pepe-Fleming weight, the censoring factor of the
# variance of the two Kaplan-Meier estimates (Pepe and Fleming, 1989).
wkm_reference <- function(time, event, grp, weight = "PF") {
  n <- length(time)
  n1 <- sum(grp == 1)
  n2 <- sum(grp == 0)
  d1 <- data.frame(t = time[grp == 1], s = event[grp == 1])
  d2 <- data.frame(t = time[grp == 0], s = event[grp == 0])

  f1 <- summary(survival::survfit(survival::Surv(t, s == 1) ~ 1, data = d1))
  f2 <- summary(survival::survfit(survival::Surv(t, s == 1) ~ 1, data = d2))
  s1_fun <- approxfun(f1$time, f1$surv, method = "constant",
                      yleft = 1, rule = 2, f = 0)
  s2_fun <- approxfun(f2$time, f2$surv, method = "constant",
                      yleft = 1, rule = 2, f = 0)

  et <- sort(time)
  grid <- (c(0, et[-n]) + et) / 2
  width <- diff(c(0, et))

  km_s1 <- s1_fun(grid)
  km_s2 <- s2_fun(grid)

  c1 <- summary(survival::survfit(survival::Surv(t, s == 0) ~ 1, data = d1))
  c2 <- summary(survival::survfit(survival::Surv(t, s == 0) ~ 1, data = d2))
  g1_fun <- if (length(c1$time) >= 1) {
    approxfun(c(c1$time, max(c1$time) + 1), c(1, c1$surv), method = "constant",
              yleft = 1, yright = min(c1$surv), rule = 2, f = 1)
  } else {
    function(x) rep(1, length(x))
  }
  g2_fun <- if (length(c2$time) >= 1) {
    approxfun(c(c2$time, max(c2$time) + 1), c(1, c2$surv), method = "constant",
              yleft = 1, yright = min(c2$surv), rule = 2, f = 1)
  } else {
    function(x) rep(1, length(x))
  }
  g1 <- g1_fun(grid)
  g2 <- g2_fun(grid)
  pf <- ifelse(g1 + g2 == 0, 0, (n * g1 * g2) / (n1 * g1 + n2 * g2))
  wt <- switch(weight, PF = pf, sqrtPF = sqrt(pf), constant = rep(1, n))

  num_raw <- sum(wt * (km_s1 - km_s2) * width)

  fp <- summary(survival::survfit(survival::Surv(time, event) ~ 1))
  sp_fun <- approxfun(fp$time, fp$surv, method = "constant",
                      yleft = 1, rule = 2, f = 0)
  sp <- sp_fun(grid)
  a_seq <- cumsum(width * wt * sp)
  a_rem <- a_seq[n] - a_seq
  sp_next <- c(sp[-1], min(fp$surv))
  d_sm <- sp_next - sp
  keep <- pf > 0 & sp > 0
  variance <- -sum((a_rem[keep]^2 / (sp[keep]^2 * pf[keep])) * d_sm[keep])

  num <- sqrt(n1 * n2 / n) * num_raw
  z <- num / sqrt(variance)
  list(num_raw = num_raw, variance = variance, z = z, p = 1 - pnorm(z))
}

test_that("wkm_fast returns the expected structure", {
  set.seed(1)
  n <- 200
  g <- rep(0:1, each = n / 2)
  tt <- c(rexp(n / 2, log(2) / 12), rexp(n / 2, log(2) / 16))
  cc <- rexp(n, rate = 0.02)
  time <- pmin(tt, cc)
  event <- as.integer(tt <= cc)

  res <- wkm_fast(time, event, group = g, control = 0)
  expect_s3_class(res, "wkm_fast")
  expect_true(all(c("wdiff", "se", "lower", "upper", "z", "chisq", "p") %in%
                    names(res)))
  expect_equal(unname(res["chisq"]), unname(res["z"])^2, tolerance = 1e-8)
  expect_true(is.finite(res["p"]) && res["p"] >= 0 && res["p"] <= 1)
})

test_that("wkm_fast matches the survfit-based reference (PF weight)", {
  skip_if_not_installed("survival")
  set.seed(11)
  n_per <- 200
  g <- rep(0:1, each = n_per)
  tt <- c(rexp(n_per, log(2) / 12), rexp(n_per, log(2) / 16))
  cc <- rexp(2 * n_per, rate = 0.02)
  time <- pmin(tt, cc)
  event <- as.integer(tt <= cc)
  fast <- wkm_fast(time, event, group = g, control = 0, side = 1,
                   weight = "PF")
  ref <- wkm_reference(time, event, g)
  expect_equal(unname(fast["z"]), ref$z, tolerance = 1e-6)
  expect_equal(unname(fast["p"]), ref$p, tolerance = 1e-6)
})

test_that("a treatment benefit gives a positive weighted difference and z", {
  set.seed(7)
  n_per <- 300
  g <- rep(0:1, each = n_per)
  tt <- c(rexp(n_per, log(2) / 10), rexp(n_per, log(2) / 18))
  cc <- rexp(2 * n_per, rate = 0.02)
  time <- pmin(tt, cc)
  event <- as.integer(tt <= cc)

  res <- wkm_fast(time, event, group = g, control = 0, side = 1)
  expect_gt(unname(res["wdiff"]), 0)
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

  res2 <- wkm_fast(time, event, group = g, control = 0, side = 2)
  res1 <- wkm_fast(time, event, group = g, control = 0, side = 1)
  z <- unname(res2["z"])
  expect_equal(unname(res1["p"]), pnorm(z, lower.tail = FALSE), tolerance = 1e-8)
  expect_equal(unname(res2["p"]), 2 * pnorm(-abs(z)), tolerance = 1e-8)
})

test_that("weight options run and give finite statistics", {
  set.seed(9)
  n_per <- 200
  g <- rep(0:1, each = n_per)
  tt <- c(rexp(n_per, log(2) / 12), rexp(n_per, log(2) / 15))
  cc <- rexp(2 * n_per, rate = 0.02)
  time <- pmin(tt, cc)
  event <- as.integer(tt <= cc)

  for (w in c("PF", "sqrtPF", "constant")) {
    res <- wkm_fast(time, event, group = g, control = 0, weight = w)
    expect_true(is.finite(res["z"]))
    expect_true(is.finite(res["wdiff"]))
  }
})

test_that("input validation works", {
  expect_error(wkm_fast(1:5, c(0, 1, 0, 1), rep(0:1, 2)), "same length")
  expect_error(
    wkm_fast(1:6, rep(2L, 6), rep(0:1, 3), control = 0), "coded as 0"
  )
  expect_error(
    wkm_fast(1:6, rep(0:1, 3), rep(1:3, 2), control = 1), "two distinct"
  )
  expect_error(
    wkm_fast(1:6, rep(0:1, 3), rep(0:1, 3)), "control must be specified"
  )
  expect_error(
    wkm_fast(1:6, rep(0:1, 3), rep(0:1, 3), control = 0, weight = "foo"),
    "should be one of"
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
    res <- wkm_fast(time, event, group = g, control = 0, side = 2)
    if (is.finite(res["p"]) && res["p"] < 0.05) reject <- reject + 1L
  }
  expect_lt(reject / nsim, 0.12)
})

test_that("wkm_fast matches the survfit-based reference for every weight", {
  skip_if_not_installed("survival")
  set.seed(12)
  n_per <- 150
  g <- rep(0:1, each = n_per)
  tt <- c(rexp(n_per, log(2) / 12), rexp(n_per, log(2) / 16))
  cc <- pmin(runif(2 * n_per, 12, 36), rexp(2 * n_per, rate = 0.02))
  time <- pmin(tt, cc)
  event <- as.integer(tt <= cc)
  for (w in c("PF", "sqrtPF", "constant")) {
    fast <- wkm_fast(time, event, group = g, control = 0, side = 1, weight = w)
    ref  <- wkm_reference(time, event, g, weight = w)
    expect_equal(unname(fast["wdiff"]), ref$num_raw, tolerance = 1e-8)
    expect_equal(unname(fast["z"]), ref$z, tolerance = 1e-6)
  }
})

test_that("wkm_fast: the standard error is calibrated for every weight", {
  # Under the null, the empirical SD of the weighted difference over repeated
  # trials should match the mean estimated SE. With the weight-specific
  # divisor of versions before 1.1.0 the ratio was about 1.15 for the constant
  # weight in this setting (independent Python simulation).
  set.seed(2026)
  nsim <- 600
  wts  <- c("PF", "sqrtPF", "constant")
  est  <- matrix(NA_real_, nsim, 3, dimnames = list(NULL, wts))
  se   <- est
  for (s in seq_len(nsim)) {
    g  <- rep(0:1, each = 150)
    tt <- rexp(300, log(2) / 12)
    cc <- pmin(runif(300, 12, 36), rexp(300, rate = 0.02))
    time  <- pmin(tt, cc)
    event <- as.integer(tt <= cc)
    for (w in wts) {
      r <- wkm_fast(time, event, group = g, control = 0, weight = w)
      est[s, w] <- r[["wdiff"]]
      se[s, w]  <- r[["se"]]
    }
  }
  ratio <- apply(est, 2, sd) / colMeans(se)
  expect_true(all(abs(ratio - 1) < 0.08), info = paste(round(ratio, 3),
                                                       collapse = ", "))
})

test_that("wkm_fast: presorted = TRUE checks the order", {
  expect_error(wkm_fast(c(3, 1, 2, 4), c(1, 1, 0, 1), c(0, 1, 0, 1),
                        control = 0, presorted = TRUE), "presorted = FALSE")
})

test_that("wkm_fast: an infinite observed time gives NA statistics", {
  time  <- c(2, 3, 5, 7, 11, 13, Inf, Inf)
  event <- c(1, 1, 1, 0, 1, 1, 0, 0)
  group <- rep(1:2, 4)
  r <- wkm_fast(time, event, group, control = 1)
  # The weighted difference comes from the C++ core without arithmetic, so
  # it is NA and not NaN.
  expect_true(is.na(r[["wdiff"]]) && !is.nan(r[["wdiff"]]))
  expect_true(all(is.na(unclass(r))))
  expect_equal(attr(r, "n"), c(4, 4))
  # A finite time in place of Inf gives a finite weighted difference.
  f <- wkm_fast(replace(time, 7:8, 20), event, group, control = 1)
  expect_true(is.finite(f[["wdiff"]]))
})

test_that("wkm_fast: time and event must be numeric and time non-negative", {
  set.seed(1)
  tt <- rexp(40, 0.1)
  ee <- rbinom(40, 1, 0.7)
  gg <- rep(0:1, 20)
  expect_error(wkm_fast(as.character(tt), ee, gg, control = 0), "numeric")
  expect_error(wkm_fast(c(-1, tt[-1]), ee, gg, control = 0), "non-negative")
  expect_error(wkm_fast(tt, factor(ee), gg, control = 0),
               "numeric or logical")
})

test_that("wkm_fast: the formula interface gives the vector result", {
  set.seed(70)
  n <- 120
  d <- data.frame(tt = round(rexp(n, 0.1), 1), ee = rbinom(n, 1, 0.8),
                  g = rep(c(1, 0), n / 2))
  ref <- wkm_fast(d$tt, d$ee, d$g, control = 0, weight = "sqrtPF")
  expect_identical(wkm_fast(Surv(tt, ee) ~ g, data = d, weight = "sqrtPF"),
                   ref)
  expect_identical(class(ref), c("wkm_fast", "fastsurv"))
  expect_error(wkm_fast(d$tt, d$ee, d$g, 0, data = d), "only with a formula")
})

test_that("wkm_fast: the formula interface passes every argument", {
  set.seed(60)
  n <- 120
  d <- data.frame(tt = round(rexp(n, 0.1), 1) + 0.1, ee = rbinom(n, 1, 0.8),
                  g = rep(c("b", "a"), n / 2), s = rep(c("x", "y", "z"), n / 3))
  d1 <- d[order(d$tt), ]
  d2 <- d[order(d$s, d$tt), ]
  expect_identical(
    wkm_fast(Surv(tt, ee) ~ g, d1, control = "b", side = 1,
             conf.level = 0.9, weight = "constant", presorted = TRUE),
    wkm_fast(d1$tt, d1$ee, d1$g, "b", side = 1, conf.level = 0.9,
             weight = "constant", presorted = TRUE))
})
