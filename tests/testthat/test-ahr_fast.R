test_that("ahr_fast returns a well-formed object", {
  set.seed(11)
  n <- 120
  obs <- c(rexp(n, 0.1), rexp(n, 0.16))
  cens <- rexp(2 * n, 0.05)
  status <- as.integer(obs <= cens)
  obs <- pmin(obs, cens)
  group <- rep(c(0, 1), each = n)

  fit <- ahr_fast(obs, status, group, control = 0, tau = 8)

  expect_s3_class(fit, "ahr_fast")
  expect_true(all(c("ahr", "log.ahr", "se.loghr", "z", "p.value", "z.loghr",
                    "p.value.loghr", "se.theta", "lower", "upper", "theta",
                    "var.theta1", "var.theta2", "tau", "n", "groups") %in%
                    names(fit)))
  # theta1 + theta2 == 1 by construction
  expect_equal(sum(fit$theta), 1, tolerance = 1e-12)
  expect_true(is.finite(fit$ahr) && fit$ahr > 0)
  expect_true(fit$var.theta1 > 0)
  expect_true(fit$var.theta2 > 0)
  expect_lt(fit$lower, fit$ahr)
  expect_lt(fit$ahr, fit$upper)
})

# Independent pure-R reference: a literal port of the Kalbfleisch-Prentice
# k = 2 estimator and its direct variance, using survival::survfit for the
# Kaplan-Meier curves and a manual Greenwood variance of log(S).
ref_ahr <- function(time, status, group, tau) {
  lev <- sort(unique(group))
  s1 <- group == lev[1]
  s2 <- group == lev[2]
  t1 <- time[s1]; d1 <- status[s1]; n1 <- sum(s1)
  t2 <- time[s2]; d2 <- status[s2]; n2 <- sum(s2)
  n <- n1 + n2
  p1 <- n1 / n
  p2 <- n2 / n

  events <- c(t1[d1 == 1], t2[d2 == 1])
  events <- events[events <= tau]
  grid <- sort(unique(c(0, events, tau)))

  km_eval <- function(tt, dd, ng) {
    sf <- survival::survfit(survival::Surv(tt, dd) ~ 1)
    s_fun <- stats::approxfun(sf$time, sf$surv, method = "constant",
                              yleft = 1, rule = 2, f = 0)
    term <- sf$n.event / (sf$n.risk * (sf$n.risk - sf$n.event))
    term[!is.finite(term)] <- 0
    gw <- cumsum(term)
    v_fun <- stats::approxfun(sf$time, ng * gw, method = "constant",
                              yleft = 0, rule = 2, f = 0)
    list(S = s_fun(grid), phi = v_fun(grid))
  }

  f1 <- km_eval(t1, d1, n1)
  f2 <- km_eval(t2, d2, n2)
  S1 <- f1$S; phi1 <- f1$phi
  S2 <- f2$S; phi2 <- f2$phi
  m <- length(grid)

  dS1 <- S1 - c(1, S1[-m])
  dS2 <- S2 - c(1, S2[-m])
  GL <- S1[m] * S2[m]
  denom <- 1 - GL

  theta1 <- -sum(S2 * dS1) / denom
  theta2 <- 1 - theta1
  ahr <- theta2 / theta1

  A12 <- GL / p1 * sum(S1 * phi1 * dS2)
  A21 <- GL / p2 * sum(S2 * phi2 * dS1)
  C1 <- GL^2 * phi1[m] / p1
  C2 <- GL^2 * phi2[m] / p2

  B122 <- 0
  for (u in 1:m) {
    inner <- sum(S1 * phi1[pmin(1:m, u)] * dS2)
    B122 <- B122 + S1[u] * inner * dS2[u]
  }
  B122 <- B122 / p1

  B211 <- 0
  for (u in 1:m) {
    inner <- sum(S2 * phi2[pmin(1:m, u)] * dS1)
    B211 <- B211 + S2[u] * inner * dS1[u]
  }
  B211 <- B211 / p2

  Vx11 <- B122 + B211 + C1 - 2 * A12
  VxG1 <- A12 - A21 - C1
  VG <- C1 + C2
  Sigma11 <- (Vx11 + 2 * theta1 * VxG1 + theta1^2 * VG) / denom^2
  var.theta1 <- Sigma11 / n

  Vx22 <- B211 + B122 + C2 - 2 * A21
  VxG2 <- A21 - A12 - C2
  Sigma22 <- (Vx22 + 2 * theta2 * VxG2 + theta2^2 * VG) / denom^2
  var.theta2 <- Sigma22 / n

  list(theta1 = theta1, theta2 = theta2, ahr = ahr,
       var.theta1 = var.theta1, var.theta2 = var.theta2)
}

test_that("point estimate matches the independent reference", {
  skip_if_not_installed("survival")
  set.seed(123)
  ne <- 80
  o1 <- rexp(ne, 0.10); c1 <- rexp(ne, 0.05)
  o2 <- rexp(ne, 0.18); c2 <- rexp(ne, 0.05)
  obs <- c(pmin(o1, c1), pmin(o2, c2))
  status <- as.integer(c(o1 <= c1, o2 <= c2))
  group <- rep(c(0, 1), each = ne)
  tau <- 8

  fit <- ahr_fast(obs, status, group, control = 0, tau = tau)
  ref <- ref_ahr(obs, status, group, tau)

  expect_equal(unname(fit$theta[1]), ref$theta1, tolerance = 1e-6)
  expect_equal(unname(fit$theta[2]), ref$theta2, tolerance = 1e-6)
  expect_equal(fit$ahr, ref$ahr, tolerance = 1e-6)
})

test_that("direct variance matches the literal reference", {
  skip_if_not_installed("survival")
  set.seed(123)
  ne <- 80
  o1 <- rexp(ne, 0.10); c1 <- rexp(ne, 0.05)
  o2 <- rexp(ne, 0.18); c2 <- rexp(ne, 0.05)
  obs <- c(pmin(o1, c1), pmin(o2, c2))
  status <- as.integer(c(o1 <= c1, o2 <= c2))
  group <- rep(c(0, 1), each = ne)
  tau <- 8

  fit <- ahr_fast(obs, status, group, control = 0, tau = tau)
  ref <- ref_ahr(obs, status, group, tau)

  expect_equal(fit$var.theta1, ref$var.theta1, tolerance = 1e-6)
  expect_equal(fit$var.theta2, ref$var.theta2, tolerance = 1e-6)
})

test_that("inference quantities are internally consistent", {
  set.seed(7)
  n <- 100
  obs <- c(rexp(n, 0.1), rexp(n, 0.15))
  cens <- rexp(2 * n, 0.04)
  status <- as.integer(obs <= cens)
  obs <- pmin(obs, cens)
  group <- rep(c(0, 1), each = n)

  fit <- ahr_fast(obs, status, group, control = 0, tau = 7)

  # primary test on the theta (comparison-share) scale, using var.theta2
  expect_equal(fit$se.theta, sqrt(fit$var.theta2), tolerance = 1e-12)
  expect_equal(fit$z, (fit$theta[[2]] - 0.5) / sqrt(fit$var.theta2),
               tolerance = 1e-10)
  expect_equal(fit$p.value, 2 * pnorm(-abs(fit$z)), tolerance = 1e-12)
  # equivalent log(AHR)-scale quantities
  expect_equal(fit$z.loghr, fit$log.ahr / fit$se.loghr, tolerance = 1e-10)
  expect_equal(fit$p.value.loghr, 2 * pnorm(-abs(fit$z.loghr)), tolerance = 1e-12)
  # theta-scale and log-scale Z agree in sign (same direction of effect)
  expect_equal(sign(fit$z), sign(fit$z.loghr))
  zq <- qnorm(1 - (1 - fit$conf.level) / 2)
  expect_equal(fit$lower, exp(fit$log.ahr - zq * fit$se.loghr), tolerance = 1e-10)
  expect_equal(fit$upper, exp(fit$log.ahr + zq * fit$se.loghr), tolerance = 1e-10)

  # null.ahr at the point estimate gives Z = 0 and p = 1 on both scales
  fit0 <- ahr_fast(obs, status, group, control = 0, tau = 7, null.ahr = fit$ahr)
  expect_equal(fit0$z, 0, tolerance = 1e-10)
  expect_equal(fit0$p.value, 1, tolerance = 1e-10)
  expect_equal(fit0$z.loghr, 0, tolerance = 1e-10)
})

test_that("estimate recovers the true hazard ratio under proportional hazards", {
  set.seed(2024)
  ne <- 3000
  rate0 <- 0.1
  psi <- 1.8
  o1 <- rexp(ne, rate0)
  o2 <- rexp(ne, rate0 * psi)
  cens <- rexp(2 * ne, 0.02)
  obs <- pmin(c(o1, o2), cens)
  status <- as.integer(c(o1, o2) <= cens)
  group <- rep(c(0, 1), each = ne)

  fit <- ahr_fast(obs, status, group, control = 0, tau = 6)
  expect_equal(fit$ahr, psi, tolerance = 0.12)
})

test_that("swapping group labels inverts the average hazard ratio", {
  skip_if_not_installed("survival")
  set.seed(55)
  ne <- 1500
  o1 <- rexp(ne, 0.1); o2 <- rexp(ne, 0.16)
  cens <- rexp(2 * ne, 0.03)
  obs <- pmin(c(o1, o2), cens)
  status <- as.integer(c(o1, o2) <= cens)
  group <- rep(c(0, 1), each = ne)

  fit <- ahr_fast(obs, status, group, control = 0, tau = 7)
  fit_sw <- ahr_fast(obs, status, 1 - group, control = 0, tau = 7)

  # the swapped call must still match the literal reference exactly
  ref_sw <- ref_ahr(obs, status, 1 - group, tau = 7)
  expect_equal(fit_sw$ahr, ref_sw$ahr, tolerance = 1e-6)

  # and the two estimates are reciprocal up to the discretization term
  expect_equal(fit$ahr * fit_sw$ahr, 1, tolerance = 0.02)
  expect_equal(fit$p.value, fit_sw$p.value, tolerance = 0.05)
})

test_that("presorted gives identical results", {
  set.seed(99)
  n <- 150
  obs <- c(rexp(n, 0.1), rexp(n, 0.2))
  cens <- rexp(2 * n, 0.05)
  status <- as.integer(obs <= cens)
  obs <- pmin(obs, cens)
  group <- rep(c(0, 1), each = n)

  ord <- order(obs)
  fit_default <- ahr_fast(obs, status, group, control = 0, tau = 6)
  fit_presort <- ahr_fast(obs[ord], status[ord], group[ord], control = 0, tau = 6,
                          presorted = TRUE)

  expect_equal(fit_default$ahr, fit_presort$ahr, tolerance = 1e-12)
  expect_equal(fit_default$var.theta1, fit_presort$var.theta1, tolerance = 1e-12)
})

test_that("factor group and default tau work", {
  set.seed(3)
  n <- 90
  obs <- c(rexp(n, 0.1), rexp(n, 0.15))
  cens <- rexp(2 * n, 0.05)
  status <- as.integer(obs <= cens)
  obs <- pmin(obs, cens)
  group <- factor(rep(c("A", "B"), each = n))

  fit <- ahr_fast(obs, status, group, control = "A")
  expect_s3_class(fit, "ahr_fast")
  expect_equal(fit$tau, min(max(obs[group == "A"]), max(obs[group == "B"])),
               tolerance = 1e-12)
  expect_identical(as.character(fit$groups), c("A", "B"))
})

test_that("input validation catches bad arguments", {
  set.seed(1)
  n <- 50
  obs <- rexp(2 * n, 0.1)
  status <- rep(c(0L, 1L), n)
  group <- rep(c(0, 1), each = n)

  # three distinct group values
  expect_error(ahr_fast(obs, status, rep(c(0, 1, 2), length.out = 2 * n)))
  # mismatched lengths
  expect_error(ahr_fast(obs[-1], status, group))
  # status not 0/1
  bad_status <- status
  bad_status[1] <- 2L
  expect_error(ahr_fast(obs, bad_status, group))
  # invalid confidence level
  expect_error(ahr_fast(obs, status, group, conf.level = 1.5))
  # non-positive null.ahr
  expect_error(ahr_fast(obs, status, group, null.ahr = -1))
  # missing values
  na_obs <- obs
  na_obs[1] <- NA
  expect_error(ahr_fast(na_obs, status, group))
  # control not one of the two group values
  expect_error(ahr_fast(obs, status, group, control = 9), "must be one of")
})

test_that("tied event times across groups follow the AHR package convention", {
  # Hand calculation (right-continuous S2 at each jump of S1, theta2 =
  # 1 - theta1, as in the AHR package) with all events and tau = 5:
  # A = 1..6 as reference gives ahr = 1 for B = (1, 2, 3, 5, 7, 8), whereas
  # B as reference gives 21/13, not 1. The asymmetry comes from the ties at
  # times 1, 2, 3, and 5.
  tt <- c(1, 2, 3, 4, 5, 6, 1, 2, 3, 5, 7, 8)
  ee <- rep(1L, 12)
  gg <- rep(c("A", "B"), each = 6)
  expect_equal(ahr_fast(tt, ee, gg, control = "A", tau = 5)$ahr, 1,
               tolerance = 1e-12)
  expect_equal(ahr_fast(tt, ee, gg, control = "B", tau = 5)$ahr, 21 / 13,
               tolerance = 1e-12)
})

test_that("without ties, exchanging the reference gives the reciprocal", {
  set.seed(515)
  tt <- c(rexp(80, 0.10), rexp(80, 0.07))
  ee <- rbinom(160, 1, 0.8)
  gg <- rep(0:1, each = 80)
  a1 <- ahr_fast(tt, ee, gg, control = 0, tau = 10)$ahr
  a2 <- ahr_fast(tt, ee, gg, control = 1, tau = 10)$ahr
  expect_equal(a1 * a2, 1, tolerance = 1e-10)
})

test_that("ahr_fast warns when a supplied tau exceeds the observed follow-up", {
  set.seed(74)
  tt <- c(rexp(50, 0.10), rexp(50, 0.07))
  ee <- rep(1L, 100)
  gg <- rep(0:1, each = 50)
  lim <- min(max(tt[gg == 0]), max(tt[gg == 1]))
  expect_warning(ahr_fast(tt, ee, gg, control = 0, tau = lim + 1), "exceeds")
  expect_warning(ahr_fast(tt, ee, gg, control = 0), NA)
})

test_that("ahr_fast: presorted = TRUE checks the order", {
  time  <- c(3, 1, 2, 4, 6, 5)
  event <- c(1, 1, 0, 1, 1, 0)
  group <- c(0, 1, 0, 1, 0, 1)
  expect_error(ahr_fast(time, event, group, control = 0, tau = 3, presorted = TRUE), "presorted = FALSE")
})

test_that("ahr_fast: event values are validated before conversion", {
  expect_error(ahr_fast(c(1, 2, 3, 4), c(1, 0.7, 0, 1), c(0, 0, 1, 1),
                        control = 0), "0 or 1")
})
