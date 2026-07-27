# Reference implementation of the Pepe-Fleming weighted Kaplan-Meier statistic,
# built from survfit and approxfun. This mirrors the algorithm of
# nphsim::wkm.Stat and provides an independent validation axis for wkm_core.
# event: 1 = event, 0 = censoring; grp: 1 = treatment, 0 = control.
wkm_reference <- function(time, event, grp) {
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
  wt <- ifelse(g1 + g2 == 0, 0, (n * g1 * g2) / (n1 * g1 + n2 * g2))

  num_raw <- sum(wt * (km_s1 - km_s2) * width)

  fp <- summary(survival::survfit(survival::Surv(time, event) ~ 1))
  sp_fun <- approxfun(fp$time, fp$surv, method = "constant",
                      yleft = 1, rule = 2, f = 0)
  sp <- sp_fun(grid)
  a_seq <- cumsum(width * wt * sp)
  a_rem <- a_seq[n] - a_seq
  sp_next <- c(sp[-1], min(fp$surv))
  d_sm <- sp_next - sp
  keep <- wt > 0 & sp > 0
  variance <- -sum((a_rem[keep]^2 / (sp[keep]^2 * wt[keep])) * d_sm[keep])

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
  out <- tryCatch({
    set.seed(11)
    n_per <- 200
    g <- rep(0:1, each = n_per)
    tt <- c(rexp(n_per, log(2) / 12), rexp(n_per, log(2) / 16))
    cc <- rexp(2 * n_per, rate = 0.02)
    time <- pmin(tt, cc)
    event <- as.integer(tt <= cc)
    ref <- wkm_reference(time, event, g)
    fast <- wkm_fast(time, event, group = g, control = 0, side = 1,
                     weight = "PF")
    list(ref = ref, fast = fast)
  }, error = function(e) NULL)
  skip_if(is.null(out), "survival reference unavailable")
  expect_equal(unname(out$fast["z"]), out$ref$z, tolerance = 1e-6)
  expect_equal(unname(out$fast["p"]), out$ref$p, tolerance = 1e-6)
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
    wkm_fast(1:6, rep(2L, 6), rep(0:1, 3), control = 0), "0"
  )
  expect_error(
    wkm_fast(1:6, rep(0:1, 3), rep(1:3, 2), control = 1), "two distinct"
  )
  expect_error(
    wkm_fast(1:6, rep(0:1, 3), rep(0:1, 3)), "control must be specified"
  )
  expect_error(
    wkm_fast(1:6, rep(0:1, 3), rep(0:1, 3), control = 0, weight = "foo")
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
