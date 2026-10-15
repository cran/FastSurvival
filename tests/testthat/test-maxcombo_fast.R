# ------------------------------------------------------------------ #
#  Reference max-combo implementation (pure R), used to verify the
#  C++ core and the wrapper. Returns the component signed Z-scores, the
#  correlation matrix, the statistic, and the p-value for a set of
#  Fleming-Harrington weights, using the left-continuous pooled KM.
# ------------------------------------------------------------------ #
ref_combo <- function(time, event, j, rho, gamma, side = 1,
                      abseps = 1e-6) {
  ord   <- order(time)
  time  <- time[ord]
  event <- as.integer(event[ord])
  j     <- j[ord]

  ut <- sort(unique(time[event == 1L]))
  nw <- length(rho)

  d  <- numeric(length(ut)); d1 <- numeric(length(ut))
  nn <- numeric(length(ut)); n1 <- numeric(length(ut))
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

  e1 <- d * n1 / nn
  v1 <- d * n1 * (nn - n1) * (nn - d) / (nn^2 * (nn - 1))
  v1[nn <= 1] <- 0
  keep <- nn > 1 & d > 0

  Umat <- numeric(nw)
  Vmat <- matrix(0, nw, nw)
  wlist <- vector("list", nw)
  for (a in seq_len(nw)) {
    wa <- s_minus^rho[a] * (1 - s_minus)^gamma[a]
    wlist[[a]] <- wa
    Umat[a] <- sum(wa[keep] * (d1[keep] - e1[keep]))
  }
  for (a in seq_len(nw)) {
    for (b in seq_len(nw)) {
      Vmat[a, b] <- sum(wlist[[a]][keep] * wlist[[b]][keep] * v1[keep])
    }
  }

  z_vec    <- Umat / sqrt(diag(Vmat))
  corr_mat <- stats::cov2cor(Vmat)

  if (side == 1L) {
    m_obs <- min(z_vec)
    lower <- rep(min(z_vec), nw); upper <- rep(Inf, nw)
  } else {
    m_obs <- max(abs(z_vec))
    lower <- rep(-m_obs, nw); upper <- rep(m_obs, nw)
  }

  if (nw == 1L) {
    joint <- stats::pnorm(upper) - stats::pnorm(lower)
  } else {
    joint <- mvtnorm::pmvnorm(lower = lower, upper = upper, corr = corr_mat,
                              algorithm = mvtnorm::GenzBretz(abseps = abseps))[1L]
  }

  list(z = z_vec, corr = corr_mat, U = Umat, V = Vmat,
       statistic = as.numeric(m_obs), p.value = 1 - as.numeric(joint))
}

test_that("component U and V match the pure-R reference", {
  skip_if_not_installed("survival")
  ov  <- survival::ovarian
  jj  <- as.integer(ov$rx != 1)
  rho <- c(0, 0, 1, 1); gamma <- c(0, 1, 0, 1)
  ref <- ref_combo(ov$futime, ov$fustat, jj, rho, gamma)
  cc  <- combo_logrank_core(ov$futime[order(ov$futime)],
                            as.integer(ov$fustat[order(ov$futime)]),
                            jj[order(ov$futime)],
                            rho, gamma)
  expect_equal(as.numeric(cc$U), as.numeric(ref$U), tolerance = 1e-8)
  expect_equal(as.numeric(cc$V), as.numeric(ref$V), tolerance = 1e-8)
})

test_that("first component FH(0,0) reproduces the ordinary log-rank Z", {
  skip_if_not_installed("survival")
  ov  <- survival::ovarian
  fit <- maxcombo_fast(ov$futime, ov$fustat, ov$rx, 1, side = 1,
                       rho = c(0, 0, 1, 1), gamma = c(0, 1, 0, 1))
  z_lr <- as.numeric(survdiff_fast(ov$futime, ov$fustat, ov$rx, 1, side = 1))
  expect_equal(as.numeric(attr(fit, "z")[1L]), z_lr, tolerance = 1e-8)
})

test_that("each component Z matches survdiff_fast with the same FH weight", {
  skip_if_not_installed("survival")
  ov  <- survival::ovarian
  rho <- c(0, 0, 1, 1); gamma <- c(0, 1, 0, 1)
  fit <- maxcombo_fast(ov$futime, ov$fustat, ov$rx, 1, side = 1,
                       rho = rho, gamma = gamma)
  z_combo <- as.numeric(attr(fit, "z"))
  for (k in seq_along(rho)) {
    z_single <- as.numeric(survdiff_fast(ov$futime, ov$fustat, ov$rx, 1,
                                         side = 1, weight = "fh",
                                         rho = rho[k], gamma = gamma[k]))
    expect_equal(z_combo[k], z_single, tolerance = 1e-8)
  }
})

test_that("statistic, correlation, and p-value match the pure-R reference", {
  skip_if_not_installed("survival")
  ov  <- survival::ovarian
  jj  <- as.integer(ov$rx != 1)
  rho <- c(0, 0, 1, 1); gamma <- c(0, 1, 0, 1)
  # Both p-values are randomized quasi-Monte-Carlo integrals whose error is of
  # the order of 1e-4 here (singular correlation), so the seed is fixed and the
  # absolute difference is compared with 2e-3.
  set.seed(1)
  ref <- ref_combo(ov$futime, ov$fustat, jj, rho, gamma, side = 1)
  fit <- maxcombo_fast(ov$futime, ov$fustat, ov$rx, 1, side = 1,
                       rho = rho, gamma = gamma)
  expect_equal(as.numeric(fit["statistic"]), ref$statistic, tolerance = 1e-8)
  expect_equal(as.numeric(attr(fit, "corr")), as.numeric(ref$corr),
               tolerance = 1e-8)
  expect_lt(abs(as.numeric(fit["p.value"]) - ref$p.value), 2e-3)
})

test_that("two-sided statistic and p-value match the pure-R reference", {
  skip_if_not_installed("survival")
  ov  <- survival::ovarian
  jj  <- as.integer(ov$rx != 1)
  rho <- c(0, 0, 1, 1); gamma <- c(0, 1, 0, 1)
  # Randomized integrals on both sides: fixed seed, absolute tolerance 2e-3.
  set.seed(1)
  ref <- ref_combo(ov$futime, ov$fustat, jj, rho, gamma, side = 2)
  fit <- maxcombo_fast(ov$futime, ov$fustat, ov$rx, 1, side = 2,
                       rho = rho, gamma = gamma)
  expect_equal(as.numeric(fit["statistic"]), ref$statistic, tolerance = 1e-8)
  expect_lt(abs(as.numeric(fit["p.value"]) - ref$p.value), 2e-3)
})

test_that("two-sided with two weights matches the pure-R reference", {
  skip_if_not_installed("survival")
  ov  <- survival::ovarian
  jj  <- as.integer(ov$rx != 1)
  rho <- c(0, 0); gamma <- c(0, 1)
  # Randomized integrals on both sides: fixed seed, absolute tolerance 2e-3.
  set.seed(1)
  ref <- ref_combo(ov$futime, ov$fustat, jj, rho, gamma, side = 2)
  fit <- maxcombo_fast(ov$futime, ov$fustat, ov$rx, 1, side = 2,
                       rho = rho, gamma = gamma)
  expect_equal(as.numeric(fit["statistic"]), ref$statistic, tolerance = 1e-8)
  expect_lt(abs(as.numeric(fit["p.value"]) - ref$p.value), 2e-3)
})

test_that("two-sided with three weights matches the pure-R reference", {
  skip_if_not_installed("survival")
  ov  <- survival::ovarian
  jj  <- as.integer(ov$rx != 1)
  rho <- c(0, 0, 1); gamma <- c(0, 1, 0)
  # Randomized integrals on both sides: fixed seed, absolute tolerance 2e-3.
  set.seed(1)
  ref <- ref_combo(ov$futime, ov$fustat, jj, rho, gamma, side = 2)
  fit <- maxcombo_fast(ov$futime, ov$fustat, ov$rx, 1, side = 2,
                       rho = rho, gamma = gamma)
  expect_equal(as.numeric(fit["statistic"]), ref$statistic, tolerance = 1e-8)
  expect_lt(abs(as.numeric(fit["p.value"]) - ref$p.value), 2e-3)
})

test_that("a single FH weight reduces to the one-sided weighted Z p-value", {
  skip_if_not_installed("survival")
  ov  <- survival::ovarian
  fit <- maxcombo_fast(ov$futime, ov$fustat, ov$rx, 1, side = 1,
                       rho = 0, gamma = 1)
  z   <- as.numeric(survdiff_fast(ov$futime, ov$fustat, ov$rx, 1, side = 1,
                                  weight = "fh", rho = 0, gamma = 1))
  # One weight: statistic is the Z itself (min over a single component). With
  # lower = z and upper = +Inf, the joint probability is 1 - pnorm(z), so the
  # p-value is 1 - (1 - pnorm(z)) = pnorm(z), the lower-tail probability that
  # is small when treatment is favoured (z negative).
  expect_equal(as.numeric(fit["statistic"]), z, tolerance = 1e-8)
  expect_equal(as.numeric(fit["p.value"]), pnorm(z), tolerance = 1e-8)
})

test_that("two-weight TVPACK path matches the pure-R reference", {
  skip_if_not_installed("survival")
  ov  <- survival::ovarian
  jj  <- as.integer(ov$rx != 1)
  rho <- c(0, 0); gamma <- c(0, 1)
  ref <- ref_combo(ov$futime, ov$fustat, jj, rho, gamma, side = 1)
  fit <- maxcombo_fast(ov$futime, ov$fustat, ov$rx, 1, side = 1,
                       rho = rho, gamma = gamma)
  expect_equal(as.numeric(fit["statistic"]), ref$statistic, tolerance = 1e-8)
  expect_equal(as.numeric(fit["p.value"]), ref$p.value, tolerance = 1e-4)
})

test_that("presorted = TRUE matches presorted = FALSE", {
  skip_if_not_installed("survival")
  ov  <- survival::ovarian
  f1  <- maxcombo_fast(ov$futime, ov$fustat, ov$rx, 1, side = 1)
  ord <- order(ov$futime)
  f2  <- maxcombo_fast(ov$futime[ord], ov$fustat[ord], ov$rx[ord], 1,
                       side = 1, presorted = TRUE)
  expect_equal(as.numeric(f1["statistic"]), as.numeric(f2["statistic"]),
               tolerance = 1e-10)
  expect_equal(as.numeric(attr(f1, "z")), as.numeric(attr(f2, "z")),
               tolerance = 1e-10)
})

test_that("character and factor group give the same result as numeric", {
  skip_if_not_installed("survival")
  ov    <- survival::ovarian
  g_num <- ov$rx
  g_chr <- paste0("arm_", ov$rx)
  g_fac <- factor(ov$rx)
  f_num <- maxcombo_fast(ov$futime, ov$fustat, g_num, control = 1, side = 1)
  f_chr <- maxcombo_fast(ov$futime, ov$fustat, g_chr, control = "arm_1", side = 1)
  f_fac <- maxcombo_fast(ov$futime, ov$fustat, g_fac, control = "1", side = 1)
  expect_equal(as.numeric(f_num["statistic"]), as.numeric(f_chr["statistic"]),
               tolerance = 1e-10)
  expect_equal(as.numeric(f_num["statistic"]), as.numeric(f_fac["statistic"]),
               tolerance = 1e-10)
})

test_that("attributes record rho, gamma, side, and n", {
  skip_if_not_installed("survival")
  ov  <- survival::ovarian
  fit <- maxcombo_fast(ov$futime, ov$fustat, ov$rx, 1, side = 2,
                       rho = c(0, 0, 1, 1), gamma = c(0, 1, 0, 1))
  expect_equal(attr(fit, "rho"), c(0, 0, 1, 1))
  expect_equal(attr(fit, "gamma"), c(0, 1, 0, 1))
  expect_equal(attr(fit, "side"), 2L)
  expect_equal(attr(fit, "n"), nrow(ov))
})

test_that("the correlation matrix is symmetric with unit diagonal", {
  skip_if_not_installed("survival")
  ov  <- survival::ovarian
  fit <- maxcombo_fast(ov$futime, ov$fustat, ov$rx, 1, side = 1)
  cm  <- attr(fit, "corr")
  expect_equal(unname(diag(cm)), rep(1, ncol(cm)), tolerance = 1e-10)
  expect_equal(unname(cm), unname(t(cm)), tolerance = 1e-10)
})

test_that("mismatched rho and gamma lengths raise an error", {
  skip_if_not_installed("survival")
  ov <- survival::ovarian
  expect_error(
    maxcombo_fast(ov$futime, ov$fustat, ov$rx, 1, rho = c(0, 1), gamma = 0),
    "same length"
  )
})

test_that("mismatched input lengths raise an error", {
  skip_if_not_installed("survival")
  ov <- survival::ovarian
  expect_error(
    maxcombo_fast(ov$futime, ov$fustat[-1L], ov$rx, 1),
    "same length"
  )
})

test_that("zero events raise an error", {
  expect_error(
    maxcombo_fast(c(1, 2, 3, 4), c(0, 0, 0, 0), c(0, 0, 1, 1), 0),
    "No events"
  )
})

# ------------------------------------------------------------------ #
#  Optional cross-check against an external max-combo implementation.
#  Interfaces differ across package versions, so this block is wrapped
#  in skip_if_not_installed and a tryCatch, and uses a loose tolerance
#  to tolerate algorithmic and Monte-Carlo differences in the p-value.
# ------------------------------------------------------------------ #
test_that("p-value is in a sensible range", {
  skip_if_not_installed("survival")
  ov  <- survival::ovarian
  fit <- maxcombo_fast(ov$futime, ov$fustat, ov$rx, 1, side = 2)
  p   <- as.numeric(fit["p.value"])
  expect_true(is.finite(p) && p >= 0 && p <= 1)
})

test_that("one-sided p-value matches simtrial::maxcombo when available", {
  skip_if_not_installed("survival")
  skip_if_not_installed("simtrial")
  ov  <- survival::ovarian
  rho <- c(0, 0, 1, 1); gamma <- c(0, 1, 0, 1)

  # simtrial::maxcombo expects a cut TTE data frame with columns tte, event,
  # treatment, and stratum; build one directly from ovarian. The treatment
  # group is rx == 2, mapped to "experimental" to match the default arm value.
  ext_p <- tryCatch({
    df <- data.frame(
      stratum   = "All",
      treatment = ifelse(ov$rx == 2, "experimental", "control"),
      tte       = ov$futime,
      event     = ov$fustat
    )
    res <- simtrial::maxcombo(df, rho = rho, gamma = gamma, return_corr = TRUE)
    as.numeric(res$p_value)
  }, error = function(e) NA_real_)

  skip_if(is.na(ext_p), "simtrial interface not compatible; cross-check skipped")

  # control = 1 means rx == 2 is the treatment group, matching the mapping
  # above, so the one-sided p-values use the same sign convention.
  fit <- maxcombo_fast(ov$futime, ov$fustat, ov$rx, control = 1, side = 1,
                       rho = rho, gamma = gamma)
  p   <- as.numeric(fit["p.value"])
  # Both p-values carry Monte-Carlo error from the QMC integration; the
  # covariance constructions differ algebraically but are mathematically
  # equal, so allow a loose tolerance.
  expect_equal(p, ext_p, tolerance = 0.01)
})

test_that("maxcombo_fast: presorted = TRUE checks the order", {
  time  <- c(3, 1, 2, 4, 6, 5)
  event <- c(1, 1, 0, 1, 1, 0)
  group <- c(0, 1, 0, 1, 0, 1)
  expect_error(maxcombo_fast(time, event, group, control = 0, presorted = TRUE), "presorted = FALSE")
})

# ------------------------------------------------------------------ #
#  Singular correlation of the Fleming-Harrington weights
# ------------------------------------------------------------------ #

# The weights of G(0,1) and G(1,0) add up to that of G(0,0) at every time
# ((1 - S) + S = 1), so Z00 = c1 Z01 + c2 Z10 exactly and the correlation
# matrix of these components is singular. The one-sided p-value
# 1 - P(Z_k >= m for all k) then reduces to a one-dimensional (three weights)
# or two-dimensional (four weights) integral of normal densities and
# distribution functions, evaluated here with stats::integrate() as a
# reference independent of mvtnorm. The formulas were checked against a
# quasi-Monte-Carlo integral with 1e7 points in Python (agreement within
# 2e-5). The components are ordered FH(0,0), FH(0,1), FH(1,0) and, with four
# weights, FH(1,1).
ref_singular_p <- function(m, corr) {
  S  <- corr[2:3, 2:3]
  cc <- solve(S, corr[1, 2:3])
  r  <- S[1, 2]
  sr <- sqrt(1 - r^2)
  # Lower limit of Z10 given Z01 = x; its kink at x0 splits the outer range.
  lower_y <- function(x) pmax(m, (m - cc[1] * x) / cc[2])
  x0  <- m * (1 - cc[2]) / cc[1]
  brk <- if (x0 > m) c(m, x0, Inf) else c(m, Inf)
  pieces <- function(f, rel) {
    sum(vapply(seq_len(length(brk) - 1L), function(k) {
      integrate(f, brk[k], brk[k + 1L], rel.tol = rel, abs.tol = 1e-12)$value
    }, numeric(1)))
  }
  if (nrow(corr) == 3L) {
    f <- function(x) {
      dnorm(x) * pnorm((lower_y(x) - r * x) / sr, lower.tail = FALSE)
    }
    joint <- pieces(f, 1e-10)
  } else {
    beta <- solve(S, corr[4, 2:3])
    s11  <- sqrt(1 - sum(beta * (S %*% beta)))
    inner <- function(x) {
      g <- function(y) {
        dnorm((y - r * x) / sr) / sr *
          pnorm((m - beta[1] * x - beta[2] * y) / s11, lower.tail = FALSE)
      }
      integrate(g, lower_y(x), Inf, rel.tol = 1e-10, abs.tol = 1e-13)$value
    }
    f <- function(x) dnorm(x) * vapply(x, inner, numeric(1))
    joint <- pieces(f, 1e-7)
  }
  list(p = 1 - joint, resid = 1 - sum(cc * (S %*% cc)))
}

test_that("one-sided p-values with singular weight sets match an independent integral", {
  skip_if_not_installed("survival")
  ov <- survival::ovarian
  # Three weights: TVPACK (deterministic).
  f3 <- maxcombo_fast(ov$futime, ov$fustat, ov$rx, control = 1, side = 1,
                      rho = c(0, 0, 1), gamma = c(0, 1, 0))
  m3 <- unname(f3["statistic"])
  r3 <- ref_singular_p(m3, attr(f3, "corr"))
  expect_lt(abs(r3$resid), 1e-8)
  expect_lt(abs(unname(f3["p.value"]) - r3$p), 2e-5)
  # Four weights: GenzBretz with the default maxpts, whose error is of the
  # order of 1e-4 for this singular problem.
  set.seed(1)
  f4 <- maxcombo_fast(ov$futime, ov$fustat, ov$rx, control = 1, side = 1)
  m4 <- unname(f4["statistic"])
  r4 <- ref_singular_p(m4, attr(f4, "corr"))
  expect_lt(abs(r4$resid), 1e-8)
  expect_lt(abs(unname(f4["p.value"]) - r4$p), 2e-3)
  # The reference is sensitive at these tolerances.
  expect_gt(abs(ref_singular_p(m3 + 0.05, attr(f3, "corr"))$p - r3$p), 2e-3)
  expect_gt(abs(ref_singular_p(m4 + 0.05, attr(f4, "corr"))$p - r4$p), 2e-3)
})

test_that("maxcombo_fast: the formula interface gives the vector result", {
  set.seed(68)
  n <- 120
  d <- data.frame(tt = round(rexp(n, 0.1), 1), ee = rbinom(n, 1, 0.8),
                  g = rep(c(1, 0), n / 2))
  # The quasi-Monte-Carlo p-value of four weights uses R's random numbers
  set.seed(1)
  f4 <- maxcombo_fast(Surv(tt, ee) ~ g, data = d, side = 1)
  set.seed(1)
  v4 <- maxcombo_fast(d$tt, d$ee, d$g, control = 0, side = 1)
  expect_identical(f4, v4)
  expect_identical(
    maxcombo_fast(Surv(tt, ee) ~ g, d, rho = c(0, 0), gamma = c(0, 1)),
    maxcombo_fast(d$tt, d$ee, d$g, 0, rho = c(0, 0), gamma = c(0, 1)))
  expect_identical(class(v4), c("maxcombo_fast", "fastsurv"))
  expect_error(maxcombo_fast(d$tt, d$ee, d$g, 0, data = d),
               "only with a formula")
})

test_that("maxcombo_fast: the formula interface passes every argument", {
  set.seed(60)
  n <- 120
  d <- data.frame(tt = round(rexp(n, 0.1), 1) + 0.1, ee = rbinom(n, 1, 0.8),
                  g = rep(c("b", "a"), n / 2), s = rep(c("x", "y", "z"), n / 3))
  d1 <- d[order(d$tt), ]
  d2 <- d[order(d$s, d$tt), ]
  # The quasi-Monte-Carlo p-value of four weights uses R's random numbers
  set.seed(1)
  mc_f <- maxcombo_fast(Surv(tt, ee) ~ g, d1, control = "b", side = 1,
                        presorted = TRUE, abseps = 1e-4, maxpts = 5000)
  set.seed(1)
  mc_v <- maxcombo_fast(d1$tt, d1$ee, d1$g, "b", side = 1, presorted = TRUE,
                        abseps = 1e-4, maxpts = 5000)
  expect_identical(mc_f, mc_v)
  expect_identical(
    maxcombo_fast(Surv(tt, ee) ~ g, d, control = "b", side = 1,
                  rho = c(0, 1), gamma = c(1, 0), abseps = 1e-7),
    maxcombo_fast(d$tt, d$ee, d$g, "b", side = 1, rho = c(0, 1),
                  gamma = c(1, 0), abseps = 1e-7))
})
