# ------------------------------------------------------------------ #
#  Reference robust modestly-weighted implementation (pure R), used to
#  verify the C++ core and the wrapper. Computes the two component
#  numerators, their variances, and the null covariance from the
#  left-continuous pooled Kaplan-Meier estimate, by an explicit at-risk
#  loop that is independent of the single-pass C++ core.
# ------------------------------------------------------------------ #
ref_rmw <- function(time, event, j, s_star) {
  ord   <- order(time)
  time  <- time[ord]
  event <- as.integer(event[ord])
  j     <- j[ord]

  ut <- sort(unique(time[event == 1L]))
  d  <- numeric(length(ut))
  d1 <- numeric(length(ut))
  nn <- numeric(length(ut))
  n1 <- numeric(length(ut))
  for (k in seq_along(ut)) {
    at_risk <- time >= ut[k]
    nn[k]   <- sum(at_risk)
    n1[k]   <- sum(at_risk & j == 1L)
    is_ev   <- event == 1L & time == ut[k]
    d[k]    <- sum(is_ev)
    d1[k]   <- sum(is_ev & j == 1L)
  }

  s       <- cumprod(1 - d / nn)
  s_minus <- c(1, s[-length(s)])
  w_mw    <- pmin(1 / s_minus, 1 / s_star)

  e1    <- d * n1 / nn
  var_d <- d * n1 * (nn - n1) * (nn - d) / (nn^2 * (nn - 1))
  var_d[nn <= 1] <- 0
  keep  <- nn > 1 & d > 0
  dev   <- d1 - e1

  U_lr <- sum(dev[keep])
  V_lr <- sum(var_d[keep])
  U_mw <- sum(w_mw[keep] * dev[keep])
  V_mw <- sum(w_mw[keep]^2 * var_d[keep])
  Cuv  <- sum(w_mw[keep] * var_d[keep])

  list(O1 = sum(d1[keep]), U_lr = U_lr, V_lr = V_lr,
       U_mw = U_mw, V_mw = V_mw, Cuv = Cuv,
       z_lr = U_lr / sqrt(V_lr), z_mw = U_mw / sqrt(V_mw),
       rho = Cuv / sqrt(V_lr * V_mw))
}

# Independent combined p-value from the component Z-scores and the null
# correlation: the bivariate normal probability is a one-dimensional integral
# of the normal density and distribution function, evaluated with
# stats::integrate() rather than mvtnorm. P(Z1 >= m, Z2 >= m) integrates
# dnorm(x) * P(Z2 >= m | Z1 = x) over x >= m, and P(|Z1| <= m, |Z2| <= m)
# integrates dnorm(x) * P(|Z2| <= m | Z1 = x) over |x| <= m. On the ovarian
# data (rho = 0.982) this agreed with the bivariate normal CDF of scipy to
# 2e-16 (Python check).
ref_rmw_p <- function(z_lr, z_mw, rho, side) {
  sr <- sqrt(1 - rho^2)
  if (side == 1L) {
    m <- min(z_lr, z_mw)
    f <- function(x) {
      dnorm(x) * pnorm((m - rho * x) / sr, lower.tail = FALSE)
    }
    joint <- integrate(f, m, Inf, rel.tol = 1e-10, abs.tol = 1e-13)$value
  } else {
    m <- max(abs(z_lr), abs(z_mw))
    f <- function(x) {
      dnorm(x) * (pnorm((m - rho * x) / sr) - pnorm((-m - rho * x) / sr))
    }
    joint <- integrate(f, -m, m, rel.tol = 1e-10, abs.tol = 1e-13)$value
  }
  1 - joint
}

test_that("rmw_core matches the pure-R reference", {
  skip_if_not_installed("survival")
  ov  <- survival::ovarian
  jj  <- as.integer(ov$rx != 1)
  ref <- ref_rmw(ov$futime, ov$fustat, jj, s_star = 0.5)

  ord  <- order(ov$futime)
  core <- rmw_core(ov$futime[ord], as.integer(ov$fustat[ord]), jj[ord], 0.5)
  expect_equal(core[1L], ref$O1,   tolerance = 1e-10)
  expect_equal(core[2L], ref$U_lr, tolerance = 1e-10)
  expect_equal(core[3L], ref$V_lr, tolerance = 1e-10)
  expect_equal(core[4L], ref$U_mw, tolerance = 1e-10)
  expect_equal(core[5L], ref$V_mw, tolerance = 1e-10)
  expect_equal(core[6L], ref$Cuv,  tolerance = 1e-10)
})

test_that("the log-rank component equals survdiff_fast log-rank Z", {
  skip_if_not_installed("survival")
  ov  <- survival::ovarian
  fit <- rmw_fast(ov$futime, ov$fustat, ov$rx, 1, side = 1, s_star = 0.5)
  z_lr_ref <- as.numeric(survdiff_fast(ov$futime, ov$fustat, ov$rx, 1, side = 1))
  expect_equal(attr(fit, "z")[["logrank"]], z_lr_ref, tolerance = 1e-10)
})

test_that("the component Z-scores and correlation match the pure-R reference", {
  skip_if_not_installed("survival")
  ov  <- survival::ovarian
  jj  <- as.integer(ov$rx != 1)
  ref <- ref_rmw(ov$futime, ov$fustat, jj, s_star = 0.5)
  fit <- rmw_fast(ov$futime, ov$fustat, ov$rx, 1, side = 1, s_star = 0.5)
  expect_equal(attr(fit, "z")[["logrank"]], ref$z_lr, tolerance = 1e-10)
  expect_equal(attr(fit, "z")[["mwlrt"]],   ref$z_mw, tolerance = 1e-10)
  expect_equal(attr(fit, "corr")[1L, 2L],   ref$rho,  tolerance = 1e-10)
})

test_that("s_star = 1 collapses the two components to the log-rank test", {
  skip_if_not_installed("survival")
  ov  <- survival::ovarian
  fit <- rmw_fast(ov$futime, ov$fustat, ov$rx, 1, side = 1, s_star = 1)
  z   <- attr(fit, "z")
  expect_equal(z[["mwlrt"]], z[["logrank"]], tolerance = 1e-10)
  expect_equal(attr(fit, "corr")[1L, 2L], 1, tolerance = 1e-8)
})

test_that("side = 1 statistic is the minimum component Z", {
  skip_if_not_installed("survival")
  ov  <- survival::ovarian
  fit <- rmw_fast(ov$futime, ov$fustat, ov$rx, 1, side = 1, s_star = 0.5)
  z   <- attr(fit, "z")
  expect_equal(as.numeric(fit)[1L], min(z[["logrank"]], z[["mwlrt"]]),
               tolerance = 1e-12)
})

test_that("side = 2 statistic is the maximum absolute component Z", {
  skip_if_not_installed("survival")
  ov  <- survival::ovarian
  fit <- rmw_fast(ov$futime, ov$fustat, ov$rx, 1, side = 2, s_star = 0.5)
  z   <- attr(fit, "z")
  expect_equal(as.numeric(fit)[1L], max(abs(z[["logrank"]]), abs(z[["mwlrt"]])),
               tolerance = 1e-12)
})

test_that("the p-value matches an independent one-dimensional integral", {
  skip_if_not_installed("survival")
  ov  <- survival::ovarian
  jj  <- as.integer(ov$rx != 1)
  # Component Z-scores and correlation from the pure-R reference, not from the
  # fitted object.
  ref <- ref_rmw(ov$futime, ov$fustat, jj, s_star = 0.5)
  for (sided in c(1L, 2L)) {
    fit   <- rmw_fast(ov$futime, ov$fustat, ov$rx, 1, side = sided, s_star = 0.5)
    p_ref <- ref_rmw_p(ref$z_lr, ref$z_mw, ref$rho, sided)
    expect_equal(as.numeric(fit)[2L], p_ref, tolerance = 1e-6)
  }
  # Python: 0.1689969050 (one-sided) and 0.3379938099 (two-sided).
  expect_equal(ref_rmw_p(ref$z_lr, ref$z_mw, ref$rho, 1L), 0.1689969050,
               tolerance = 1e-8)
  expect_equal(ref_rmw_p(ref$z_lr, ref$z_mw, ref$rho, 2L), 0.3379938099,
               tolerance = 1e-8)
})

test_that("the p-value lies in the unit interval", {
  skip_if_not_installed("survival")
  ov <- survival::ovarian
  p1 <- as.numeric(rmw_fast(ov$futime, ov$fustat, ov$rx, 1, side = 1,
                            s_star = 0.5))[2L]
  p2 <- as.numeric(rmw_fast(ov$futime, ov$fustat, ov$rx, 1, side = 2,
                            s_star = 0.5))[2L]
  expect_true(p1 >= 0 && p1 <= 1)
  expect_true(p2 >= 0 && p2 <= 1)
})

test_that("presorted = TRUE matches presorted = FALSE", {
  skip_if_not_installed("survival")
  ov  <- survival::ovarian
  f1  <- rmw_fast(ov$futime, ov$fustat, ov$rx, 1, side = 1, s_star = 0.5)
  ord <- order(ov$futime)
  f2  <- rmw_fast(ov$futime[ord], ov$fustat[ord], ov$rx[ord], 1, side = 1,
                  s_star = 0.5, presorted = TRUE)
  expect_equal(as.numeric(f1), as.numeric(f2), tolerance = 1e-10)
  expect_equal(attr(f1, "z"), attr(f2, "z"), tolerance = 1e-10)
})

test_that("invalid inputs raise informative errors", {
  skip_if_not_installed("survival")
  ov <- survival::ovarian
  expect_error(
    rmw_fast(ov$futime, ov$fustat, c(1, 2), 1, side = 1),
    "same length")
  expect_error(
    rmw_fast(ov$futime, ov$fustat, ov$rx, 1, side = 3),
    "side")
  expect_error(
    rmw_fast(ov$futime, ov$fustat, ov$rx, 1, side = 1, s_star = 0),
    "s_star")
  expect_error(
    rmw_fast(ov$futime, ov$fustat, ov$rx, 1, side = 1, s_star = 1.5),
    "s_star")
  expect_error(
    rmw_fast(ov$futime, rep(0, nrow(ov)), ov$rx, 1, side = 1),
    "No events")
})

test_that("the returned object carries the expected attributes", {
  skip_if_not_installed("survival")
  ov  <- survival::ovarian
  fit <- rmw_fast(ov$futime, ov$fustat, ov$rx, 1, side = 1, s_star = 0.5)
  expect_s3_class(fit, "rmw_fast")
  expect_equal(attr(fit, "s_star"), 0.5)
  expect_equal(attr(fit, "n"), nrow(ov))
  expect_equal(dim(attr(fit, "corr")), c(2L, 2L))
  expect_named(attr(fit, "z"), c("logrank", "mwlrt"))
})

test_that("zero component variance yields an NA result", {
  # Control subject is censored early, so at every event time only the
  # treatment group is at risk and both component variances are zero.
  fit <- rmw_fast(c(1, 5, 6, 7), c(0, 1, 1, 1), c(0, 1, 1, 1),
                  control = 0, side = 1, s_star = 0.5)
  expect_s3_class(fit, "rmw_fast")
  expect_true(is.na(as.numeric(fit)[1L]))
  expect_true(is.na(as.numeric(fit)[2L]))
})

test_that("rmw_fast: presorted = TRUE checks the order", {
  time  <- c(3, 1, 2, 4, 6, 5)
  event <- c(1, 1, 0, 1, 1, 0)
  group <- c(0, 1, 0, 1, 0, 1)
  expect_error(rmw_fast(time, event, group, control = 0, presorted = TRUE), "presorted = FALSE")
})

test_that("rmw_fast: the formula interface gives the vector result", {
  set.seed(69)
  n <- 120
  d <- data.frame(tt = round(rexp(n, 0.1), 1), ee = rbinom(n, 1, 0.8),
                  g = rep(c(1, 0), n / 2))
  ref <- rmw_fast(d$tt, d$ee, d$g, control = 0, side = 1, s_star = 0.6)
  expect_identical(rmw_fast(Surv(tt, ee) ~ g, data = d, side = 1,
                            s_star = 0.6), ref)
  expect_identical(class(ref), c("rmw_fast", "fastsurv"))
  expect_error(rmw_fast(d$tt, d$ee, d$g, 0, data = d), "only with a formula")
})

test_that("rmw_fast: the formula interface passes every argument", {
  set.seed(60)
  n <- 120
  d <- data.frame(tt = round(rexp(n, 0.1), 1) + 0.1, ee = rbinom(n, 1, 0.8),
                  g = rep(c("b", "a"), n / 2), s = rep(c("x", "y", "z"), n / 3))
  d1 <- d[order(d$tt), ]
  d2 <- d[order(d$s, d$tt), ]
  expect_identical(
    rmw_fast(Surv(tt, ee) ~ g, d1, control = "b", side = 1, s_star = 0.7,
             presorted = TRUE),
    rmw_fast(d1$tt, d1$ee, d1$g, "b", side = 1, s_star = 0.7,
             presorted = TRUE))
})
