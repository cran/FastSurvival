# ------------------------------------------------------------------ #
#  Reference single-group AHSW implementation (pure R), mirroring the
#  survAH::ah1 algorithm, used to verify the C++ core. Returns the
#  average hazard and the two variance sums v_Q (log scale) and v_U
#  (identity scale) for a single group.
# ------------------------------------------------------------------ #
ref_ahsw1 <- function(time, status, tau) {
  n    <- length(time)
  indx <- order(time)
  x    <- time[indx]
  d    <- as.integer(status[indx])

  timegrd <- unique(sort(c(0, x[d == 1 & x <= tau], tau)))
  k       <- length(timegrd)

  dDt  <- numeric(k)
  Ybar <- numeric(k)
  for (jj in seq_len(k)) {
    tg       <- timegrd[jj]
    dDt[jj]  <- sum(x == tg & d == 1)
    Ybar[jj] <- sum(x >= tg)
  }

  tmp        <- dDt / Ybar
  tmp[Ybar == 0] <- 0
  dHt        <- tmp                       # Nelson-Aalen increment at each grid pt
  St         <- cumprod(1 - tmp)          # right-continuous KM at grid points

  # Running RMST(t) = integral_0^t S, using left survival on each interval
  RMSTt <- numeric(k)
  for (jj in 2:k) {
    RMSTt[jj] <- as.numeric(diff(timegrd[1:jj]) %*% St[1:(jj - 1)])
  }

  Gt       <- Ybar / n
  F_tau    <- 1 - St[k]
  RMST_tau <- RMSTt[k]

  kerQ <- 1 / F_tau - RMSTt / RMST_tau
  kerU <- 1 / RMST_tau - F_tau * RMSTt / RMST_tau^2
  v_Q  <- sum((dHt * kerQ^2 / Gt)[dHt != 0])
  v_U  <- sum((dHt * kerU^2 / Gt)[dHt != 0])

  list(F_tau = F_tau, RMST_tau = RMST_tau, AH = F_tau / RMST_tau,
       v_Q = v_Q, v_U = v_U, n = n, surv_tau = St[k])
}

test_that("core average hazard and variances match the pure-R reference", {
  skip_if_not_installed("survival")
  ov  <- survival::ovarian
  tau <- 600
  for (rx in 1:2) {
    sub <- ov[ov$rx == rx, ]
    ref <- ref_ahsw1(sub$futime, sub$fustat, tau)
    ord <- order(sub$futime)
    cc  <- ahsw_core(sub$futime[ord], as.integer(sub$fustat[ord]), tau)
    expect_equal(cc[1L], ref$F_tau,    tolerance = 1e-8)   # F(tau)
    expect_equal(cc[2L], ref$RMST_tau, tolerance = 1e-8)   # RMST(tau)
    expect_equal(cc[3L], ref$AH,       tolerance = 1e-8)   # AH
    expect_equal(cc[4L], ref$v_Q,      tolerance = 1e-8)   # v_Q
    expect_equal(cc[5L], ref$v_U,      tolerance = 1e-8)   # v_U
  }
})

test_that("per-group average hazard matches the pure-R reference", {
  skip_if_not_installed("survival")
  ov  <- survival::ovarian
  tau <- 600
  fit <- ahsw_fast(ov$futime, ov$fustat, ov$rx, control = 1, tau = tau)
  ref0 <- ref_ahsw1(ov$futime[ov$rx == 1], ov$fustat[ov$rx == 1], tau)
  ref1 <- ref_ahsw1(ov$futime[ov$rx == 2], ov$fustat[ov$rx == 2], tau)
  expect_equal(as.numeric(fit["ah.ctrl"]), ref0$AH, tolerance = 1e-8)
  expect_equal(as.numeric(fit["ah.trt"]),  ref1$AH, tolerance = 1e-8)
})

test_that("RAH and DAH match the pure-R reference contrasts", {
  skip_if_not_installed("survival")
  ov  <- survival::ovarian
  tau <- 600
  fit <- ahsw_fast(ov$futime, ov$fustat, ov$rx, control = 1, tau = tau)
  ref0 <- ref_ahsw1(ov$futime[ov$rx == 1], ov$fustat[ov$rx == 1], tau)
  ref1 <- ref_ahsw1(ov$futime[ov$rx == 2], ov$fustat[ov$rx == 2], tau)

  rah_ref <- ref1$AH / ref0$AH
  dah_ref <- ref1$AH - ref0$AH
  expect_equal(as.numeric(fit["rah"]), rah_ref, tolerance = 1e-8)
  expect_equal(as.numeric(fit["dah"]), dah_ref, tolerance = 1e-8)

  # p-values from the pure-R reference
  se_rah <- sqrt(ref1$v_Q / ref1$n + ref0$v_Q / ref0$n)
  se_dah <- sqrt(ref1$v_U / ref1$n + ref0$v_U / ref0$n)
  p_rah  <- 2 * pnorm(-abs(log(rah_ref)) / se_rah)
  p_dah  <- 2 * pnorm(-abs(dah_ref) / se_dah)
  expect_equal(as.numeric(fit["p.rah"]), p_rah, tolerance = 1e-8)
  expect_equal(as.numeric(fit["p.dah"]), p_dah, tolerance = 1e-8)
})

test_that("confidence interval endpoints match the pure-R reference", {
  skip_if_not_installed("survival")
  ov  <- survival::ovarian
  tau <- 600
  conf <- 0.95
  z   <- qnorm(1 - (1 - conf) / 2)
  fit <- ahsw_fast(ov$futime, ov$fustat, ov$rx, control = 1, tau = tau,
                   conf.level = conf)
  ref0 <- ref_ahsw1(ov$futime[ov$rx == 1], ov$fustat[ov$rx == 1], tau)
  ref1 <- ref_ahsw1(ov$futime[ov$rx == 2], ov$fustat[ov$rx == 2], tau)

  log_rah <- log(ref1$AH / ref0$AH)
  se_rah  <- sqrt(ref1$v_Q / ref1$n + ref0$v_Q / ref0$n)
  expect_equal(as.numeric(fit["rah.lower"]), exp(log_rah - z * se_rah),
               tolerance = 1e-8)
  expect_equal(as.numeric(fit["rah.upper"]), exp(log_rah + z * se_rah),
               tolerance = 1e-8)

  dah    <- ref1$AH - ref0$AH
  se_dah <- sqrt(ref1$v_U / ref1$n + ref0$v_U / ref0$n)
  expect_equal(as.numeric(fit["dah.lower"]), dah - z * se_dah, tolerance = 1e-8)
  expect_equal(as.numeric(fit["dah.upper"]), dah + z * se_dah, tolerance = 1e-8)
})

test_that("presorted = TRUE matches presorted = FALSE", {
  skip_if_not_installed("survival")
  ov  <- survival::ovarian
  tau <- 600
  f1  <- ahsw_fast(ov$futime, ov$fustat, ov$rx, control = 1, tau = tau)
  ord <- order(ov$futime)
  f2  <- ahsw_fast(ov$futime[ord], ov$fustat[ord], ov$rx[ord], control = 1,
                   tau = tau, presorted = TRUE)
  expect_equal(as.numeric(f1), as.numeric(f2), tolerance = 1e-10)
})

test_that("character and factor group give the same result as numeric", {
  skip_if_not_installed("survival")
  ov    <- survival::ovarian
  tau   <- 600
  f_num <- ahsw_fast(ov$futime, ov$fustat, ov$rx, control = 1, tau = tau)
  f_chr <- ahsw_fast(ov$futime, ov$fustat, paste0("arm_", ov$rx),
                     control = "arm_1", tau = tau)
  f_fac <- ahsw_fast(ov$futime, ov$fustat, factor(ov$rx), control = "1",
                     tau = tau)
  expect_equal(as.numeric(f_num["rah"]), as.numeric(f_chr["rah"]),
               tolerance = 1e-10)
  expect_equal(as.numeric(f_num["rah"]), as.numeric(f_fac["rah"]),
               tolerance = 1e-10)
})

test_that("attributes record tau, conf.level, and control", {
  skip_if_not_installed("survival")
  ov  <- survival::ovarian
  fit <- ahsw_fast(ov$futime, ov$fustat, ov$rx, control = 1, tau = 600,
                   conf.level = 0.9)
  expect_equal(attr(fit, "tau"), 600)
  expect_equal(attr(fit, "conf.level"), 0.9)
  expect_equal(attr(fit, "control"), 1)
})

test_that("non-positive tau raises an error", {
  skip_if_not_installed("survival")
  ov <- survival::ovarian
  expect_error(
    ahsw_fast(ov$futime, ov$fustat, ov$rx, control = 1, tau = 0),
    "positive"
  )
})

test_that("mismatched input lengths raise an error", {
  skip_if_not_installed("survival")
  ov <- survival::ovarian
  expect_error(
    ahsw_fast(ov$futime, ov$fustat[-1L], ov$rx, control = 1, tau = 600),
    "same length"
  )
})

# ------------------------------------------------------------------ #
#  Cross-check against survAH::ah2 when available. The estimates and
#  p-values should agree to machine precision because the algorithm is
#  the same; a loose tolerance is used only to absorb minor numerical
#  differences in the survival/RMST construction.
# ------------------------------------------------------------------ #
test_that("estimates and p-values match survAH::ah2 when available", {
  skip_if_not_installed("survival")
  skip_if_not_installed("survAH")
  ov  <- survival::ovarian
  tau <- 600

  ext <- tryCatch({
    arm <- as.numeric(ov$rx == 2)
    survAH::ah2(time = ov$futime, status = ov$fustat, arm = arm, tau = tau,
                conf.int = 0.95)
  }, error = function(e) NULL)

  skip_if(is.null(ext), "survAH interface not compatible; cross-check skipped")

  fit <- ahsw_fast(ov$futime, ov$fustat, ov$rx, control = 1, tau = tau)

  ah0_ext <- ext$ah["AH (arm0)", "Est."]
  ah1_ext <- ext$ah["AH (arm1)", "Est."]
  rah_ext <- ext$rah[1L, "Est."]
  dah_ext <- ext$dah[1L, "Est."]
  prah_ext <- ext$rah[1L, "P-value"]
  pdah_ext <- ext$dah[1L, "P-value"]

  expect_equal(as.numeric(fit["ah.ctrl"]), as.numeric(ah0_ext), tolerance = 1e-6)
  expect_equal(as.numeric(fit["ah.trt"]),  as.numeric(ah1_ext), tolerance = 1e-6)
  expect_equal(as.numeric(fit["rah"]),     as.numeric(rah_ext), tolerance = 1e-6)
  expect_equal(as.numeric(fit["dah"]),     as.numeric(dah_ext), tolerance = 1e-6)
  expect_equal(as.numeric(fit["p.rah"]),   as.numeric(prah_ext), tolerance = 1e-6)
  expect_equal(as.numeric(fit["p.dah"]),   as.numeric(pdah_ext), tolerance = 1e-6)
})

test_that("ahsw_fast warns when tau exceeds the observed follow-up", {
  set.seed(72)
  tt <- c(rexp(50, 0.10), rexp(50, 0.07))
  ee <- rep(1L, 100)
  gg <- rep(0:1, each = 50)
  lim <- min(max(tt[gg == 0]), max(tt[gg == 1]))
  expect_warning(ahsw_fast(tt, ee, gg, control = 0, tau = lim + 1), "exceeds")
  expect_warning(ahsw_fast(tt, ee, gg, control = 0, tau = lim / 2), NA)
})

test_that("ahsw_fast: presorted = TRUE checks the order", {
  time  <- c(3, 1, 2, 4, 6, 5)
  event <- c(1, 1, 0, 1, 1, 0)
  group <- c(0, 1, 0, 1, 0, 1)
  expect_error(ahsw_fast(time, event, group, control = 0, tau = 3, presorted = TRUE), "presorted = FALSE")
})
