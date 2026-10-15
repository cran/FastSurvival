test_that("unstratified two-sided statistic matches survdiff", {
  skip_if_not_installed("survival")
  ov  <- survival::ovarian
  fit <- survdiff_fast(ov$futime, ov$fustat, ov$rx, control = 1, side = 2)
  ref <- survival::survdiff(survival::Surv(futime, fustat) ~ rx, data = ov)
  expect_equal(as.numeric(fit), ref$chisq, tolerance = 1e-8)
})

test_that("strata = NULL is identical to omitting strata", {
  skip_if_not_installed("survival")
  ov <- survival::ovarian
  f_default <- survdiff_fast(ov$futime, ov$fustat, ov$rx, control = 1, side = 2)
  f_null    <- survdiff_fast(ov$futime, ov$fustat, ov$rx, control = 1, side = 2,
                             strata = NULL)
  expect_identical(f_default, f_null)
})

test_that("unstratified one-sided Z squares to the two-sided statistic", {
  skip_if_not_installed("survival")
  ov  <- survival::ovarian
  z   <- as.numeric(survdiff_fast(ov$futime, ov$fustat, ov$rx, 1, side = 1))
  chi <- as.numeric(survdiff_fast(ov$futime, ov$fustat, ov$rx, 1, side = 2))
  expect_equal(z^2, chi, tolerance = 1e-10)
})

test_that("stratified chi-square matches survdiff with strata()", {
  skip_if_not_installed("survival")
  ov  <- survival::ovarian
  fit <- survdiff_fast(ov$futime, ov$fustat, ov$rx, control = 1, side = 2,
                       strata = ov$resid.ds)
  ref <- survival::survdiff(
    survival::Surv(futime, fustat) ~ rx + survival::strata(resid.ds),
    data = ov)
  expect_equal(as.numeric(fit), ref$chisq, tolerance = 1e-8)
})

test_that("stratified observed and expected match survdiff (treatment group)", {
  skip_if_not_installed("survival")
  ov  <- survival::ovarian
  fit <- survdiff_fast(ov$futime, ov$fustat, ov$rx, control = 1, side = 1,
                       strata = ov$resid.ds)
  ref <- survival::survdiff(
    survival::Surv(futime, fustat) ~ rx + survival::strata(resid.ds),
    data = ov)
  # With a strata() term survdiff returns obs and exp as group-by-stratum
  # matrices, so sum over strata to obtain the per-group totals. Group levels
  # sort to (1 = control, 2 = treatment), so the treatment group is the second
  # element of the summed vectors.
  obs_grp <- if (is.matrix(ref$obs)) rowSums(ref$obs) else ref$obs
  exp_grp <- if (is.matrix(ref$exp)) rowSums(ref$exp) else ref$exp
  expect_equal(attr(fit, "O1"), as.numeric(obs_grp[2L]), tolerance = 1e-8)
  expect_equal(attr(fit, "E1"), as.numeric(exp_grp[2L]), tolerance = 1e-8)
})

test_that("a single stratum reduces to the unstratified test", {
  skip_if_not_installed("survival")
  ov  <- survival::ovarian
  one <- rep(1L, nrow(ov))
  f_str   <- survdiff_fast(ov$futime, ov$fustat, ov$rx, 1, side = 2, strata = one)
  f_plain <- survdiff_fast(ov$futime, ov$fustat, ov$rx, 1, side = 2)
  expect_equal(as.numeric(f_str), as.numeric(f_plain), tolerance = 1e-10)
})

test_that("stratified presorted = TRUE matches presorted = FALSE", {
  skip_if_not_installed("survival")
  ov  <- survival::ovarian
  f1  <- survdiff_fast(ov$futime, ov$fustat, ov$rx, 1, side = 2,
                       strata = ov$resid.ds)
  ord <- order(ov$resid.ds, ov$futime)
  f2  <- survdiff_fast(ov$futime[ord], ov$fustat[ord], ov$rx[ord], 1, side = 2,
                       presorted = TRUE, strata = ov$resid.ds[ord])
  expect_equal(as.numeric(f1), as.numeric(f2), tolerance = 1e-10)
})

test_that("stratified one-sided Z squares to the two-sided statistic", {
  skip_if_not_installed("survival")
  ov  <- survival::ovarian
  z   <- as.numeric(survdiff_fast(ov$futime, ov$fustat, ov$rx, 1, side = 1,
                                  strata = ov$resid.ds))
  chi <- as.numeric(survdiff_fast(ov$futime, ov$fustat, ov$rx, 1, side = 2,
                                  strata = ov$resid.ds))
  expect_equal(z^2, chi, tolerance = 1e-10)
})

test_that("the strata attribute records the number of strata", {
  skip_if_not_installed("survival")
  ov  <- survival::ovarian
  fit <- survdiff_fast(ov$futime, ov$fustat, ov$rx, 1, side = 2,
                       strata = ov$resid.ds)
  expect_equal(attr(fit, "strata"), length(unique(ov$resid.ds)))
})

test_that("strata of mismatched length raises an error", {
  skip_if_not_installed("survival")
  ov <- survival::ovarian
  expect_error(
    survdiff_fast(ov$futime, ov$fustat, ov$rx, 1, side = 2,
                  strata = c(1L, 2L)),
    "same length"
  )
})

test_that("character and factor strata give the same result as integer", {
  skip_if_not_installed("survival")
  ov   <- survival::ovarian
  s_int <- ov$resid.ds
  s_chr <- paste0("stratum_", ov$resid.ds)
  s_fac <- factor(ov$resid.ds)
  f_int <- survdiff_fast(ov$futime, ov$fustat, ov$rx, 1, side = 2, strata = s_int)
  f_chr <- survdiff_fast(ov$futime, ov$fustat, ov$rx, 1, side = 2, strata = s_chr)
  f_fac <- survdiff_fast(ov$futime, ov$fustat, ov$rx, 1, side = 2, strata = s_fac)
  expect_equal(as.numeric(f_int), as.numeric(f_chr), tolerance = 1e-10)
  expect_equal(as.numeric(f_int), as.numeric(f_fac), tolerance = 1e-10)
})

# ------------------------------------------------------------------ #
#  Reference weighted log-rank implementation (pure R), used to verify
#  the C++ core. Computes the signed Z = U / sqrt(V) for the treatment
#  group, with weights derived from the left-continuous pooled KM.
# ------------------------------------------------------------------ #
ref_wlr_z <- function(time, event, j, scheme, rho = 0, gamma = 0, t_star = 0) {
  ord   <- order(time)
  time  <- time[ord]
  event <- as.integer(event[ord])
  j     <- j[ord]

  ut <- sort(unique(time[event == 1L]))
  n  <- length(time)

  # Pooled risk table at distinct event times
  d  <- numeric(length(ut))
  d1 <- numeric(length(ut))
  nn <- numeric(length(ut))
  n1 <- numeric(length(ut))
  for (k in seq_along(ut)) {
    at_risk    <- time >= ut[k]
    nn[k]      <- sum(at_risk)
    n1[k]      <- sum(at_risk & j == 1L)
    is_ev      <- event == 1L & time == ut[k]
    d[k]       <- sum(is_ev)
    d1[k]      <- sum(is_ev & j == 1L)
  }

  # Left-continuous pooled KM
  s       <- cumprod(1 - d / nn)
  s_minus <- c(1, s[-length(s)])

  w <- switch(as.character(scheme),
    "0" = s_minus^rho * (1 - s_minus)^gamma,
    "1" = {
      if (t_star <= 0) {
        rep(1, length(ut))
      } else {
        pre <- s[ut < t_star]
        mw  <- if (length(pre) > 0) 1 / pre[length(pre)] else 1
        pmin(1 / s_minus, mw)
      }
    },
    "2" = nn,
    "3" = sqrt(nn)
  )

  e1 <- d * n1 / nn
  v1 <- d * n1 * (nn - n1) * (nn - d) / (nn^2 * (nn - 1))
  v1[nn <= 1] <- 0

  keep <- nn > 1 & d > 0
  U <- sum(w[keep] * (d1[keep] - e1[keep]))
  V <- sum(w[keep]^2 * v1[keep])
  U / sqrt(V)
}

test_that("weight = 'logrank' is identical to the default", {
  skip_if_not_installed("survival")
  ov <- survival::ovarian
  f_default <- survdiff_fast(ov$futime, ov$fustat, ov$rx, 1, side = 2)
  f_lr      <- survdiff_fast(ov$futime, ov$fustat, ov$rx, 1, side = 2,
                             weight = "logrank")
  expect_identical(f_default, f_lr)
})

test_that("Fleming-Harrington G(0,0) reduces to the ordinary log-rank", {
  skip_if_not_installed("survival")
  ov  <- survival::ovarian
  f_lr <- survdiff_fast(ov$futime, ov$fustat, ov$rx, 1, side = 2)
  f_fh <- survdiff_fast(ov$futime, ov$fustat, ov$rx, 1, side = 2,
                        weight = "fh", rho = 0, gamma = 0)
  expect_equal(as.numeric(f_fh), as.numeric(f_lr), tolerance = 1e-10)
})

test_that("Fleming-Harrington G(rho,0) matches survdiff with rho", {
  skip_if_not_installed("survival")
  ov  <- survival::ovarian
  f_fh <- survdiff_fast(ov$futime, ov$fustat, ov$rx, 1, side = 2,
                        weight = "fh", rho = 1, gamma = 0)
  ref  <- survival::survdiff(survival::Surv(futime, fustat) ~ rx,
                             data = ov, rho = 1)
  expect_equal(as.numeric(f_fh), ref$chisq, tolerance = 1e-6)
})

test_that("Fleming-Harrington G(0,1) matches the pure-R reference", {
  skip_if_not_installed("survival")
  ov  <- survival::ovarian
  jj  <- as.integer(ov$rx != 1)
  z_ref <- ref_wlr_z(ov$futime, ov$fustat, jj, scheme = 0, rho = 0, gamma = 1)
  z     <- as.numeric(survdiff_fast(ov$futime, ov$fustat, ov$rx, 1, side = 1,
                                    weight = "fh", rho = 0, gamma = 1))
  expect_equal(z, z_ref, tolerance = 1e-8)
})

test_that("modestly-weighted test matches the pure-R reference", {
  skip_if_not_installed("survival")
  ov  <- survival::ovarian
  jj  <- as.integer(ov$rx != 1)
  z_ref <- ref_wlr_z(ov$futime, ov$fustat, jj, scheme = 1, t_star = 365)
  z     <- as.numeric(survdiff_fast(ov$futime, ov$fustat, ov$rx, 1, side = 1,
                                    weight = "mwlrt", t_star = 365))
  expect_equal(z, z_ref, tolerance = 1e-8)
})

test_that("modestly-weighted test with t_star = 0 reduces to the log-rank", {
  skip_if_not_installed("survival")
  ov   <- survival::ovarian
  z_lr <- as.numeric(survdiff_fast(ov$futime, ov$fustat, ov$rx, 1, side = 1))
  z_mw <- as.numeric(survdiff_fast(ov$futime, ov$fustat, ov$rx, 1, side = 1,
                                   weight = "mwlrt", t_star = 0))
  expect_equal(z_mw, z_lr, tolerance = 1e-10)
})

test_that("Gehan-Breslow and Tarone-Ware match the pure-R reference", {
  skip_if_not_installed("survival")
  ov  <- survival::ovarian
  jj  <- as.integer(ov$rx != 1)
  z_gehan_ref <- ref_wlr_z(ov$futime, ov$fustat, jj, scheme = 2)
  z_tw_ref    <- ref_wlr_z(ov$futime, ov$fustat, jj, scheme = 3)
  z_gehan <- as.numeric(survdiff_fast(ov$futime, ov$fustat, ov$rx, 1, side = 1,
                                      weight = "gehan"))
  z_tw    <- as.numeric(survdiff_fast(ov$futime, ov$fustat, ov$rx, 1, side = 1,
                                      weight = "tarone-ware"))
  expect_equal(z_gehan, z_gehan_ref, tolerance = 1e-8)
  expect_equal(z_tw,    z_tw_ref,    tolerance = 1e-8)
})

test_that("weighted presorted = TRUE matches presorted = FALSE", {
  skip_if_not_installed("survival")
  ov  <- survival::ovarian
  f1  <- survdiff_fast(ov$futime, ov$fustat, ov$rx, 1, side = 1,
                       weight = "fh", rho = 0, gamma = 1)
  ord <- order(ov$futime)
  f2  <- survdiff_fast(ov$futime[ord], ov$fustat[ord], ov$rx[ord], 1, side = 1,
                       presorted = TRUE, weight = "fh", rho = 0, gamma = 1)
  expect_equal(as.numeric(f1), as.numeric(f2), tolerance = 1e-10)
})

test_that("weighted one-sided Z squares to the two-sided statistic", {
  skip_if_not_installed("survival")
  ov  <- survival::ovarian
  z   <- as.numeric(survdiff_fast(ov$futime, ov$fustat, ov$rx, 1, side = 1,
                                  weight = "fh", rho = 0, gamma = 1))
  chi <- as.numeric(survdiff_fast(ov$futime, ov$fustat, ov$rx, 1, side = 2,
                                  weight = "fh", rho = 0, gamma = 1))
  expect_equal(z^2, chi, tolerance = 1e-10)
})

test_that("mwlrt without t_star raises an error", {
  skip_if_not_installed("survival")
  ov <- survival::ovarian
  expect_error(
    survdiff_fast(ov$futime, ov$fustat, ov$rx, 1, side = 1, weight = "mwlrt"),
    "t_star"
  )
})

test_that("the weight attribute records the scheme name", {
  skip_if_not_installed("survival")
  ov  <- survival::ovarian
  fit <- survdiff_fast(ov$futime, ov$fustat, ov$rx, 1, side = 1,
                       weight = "fh", rho = 0, gamma = 1)
  expect_equal(attr(fit, "weight"), "fh")
})

# ------------------------------------------------------------------ #
#  Reference stratified weighted log-rank: per-stratum U and V from the
#  pure-R weighted reference, summed across strata and standardized once.
# ------------------------------------------------------------------ #
ref_swlr_z <- function(time, event, j, strata,
                       scheme, rho = 0, gamma = 0, t_star = 0) {
  uv_one <- function(idx) {
    tt <- time[idx]; ee <- event[idx]; jj <- j[idx]
    ord <- order(tt)
    tt <- tt[ord]; ee <- as.integer(ee[ord]); jj <- jj[ord]

    ut <- sort(unique(tt[ee == 1L]))
    if (length(ut) == 0L) return(c(U = 0, V = 0))

    d  <- numeric(length(ut)); d1 <- numeric(length(ut))
    nn <- numeric(length(ut)); n1 <- numeric(length(ut))
    for (k in seq_along(ut)) {
      at_risk <- tt >= ut[k]
      nn[k]   <- sum(at_risk)
      n1[k]   <- sum(at_risk & jj == 1L)
      is_ev   <- ee == 1L & tt == ut[k]
      d[k]    <- sum(is_ev)
      d1[k]   <- sum(is_ev & jj == 1L)
    }

    s       <- cumprod(1 - d / nn)
    s_minus <- c(1, s[-length(s)])

    w <- switch(as.character(scheme),
      "0" = s_minus^rho * (1 - s_minus)^gamma,
      "1" = {
        if (t_star <= 0) {
          rep(1, length(ut))
        } else {
          pre <- s[ut < t_star]
          mw  <- if (length(pre) > 0) 1 / pre[length(pre)] else 1
          pmin(1 / s_minus, mw)
        }
      },
      "2" = nn,
      "3" = sqrt(nn)
    )

    e1 <- d * n1 / nn
    v1 <- d * n1 * (nn - n1) * (nn - d) / (nn^2 * (nn - 1))
    v1[nn <= 1] <- 0
    keep <- nn > 1 & d > 0
    c(U = sum(w[keep] * (d1[keep] - e1[keep])),
      V = sum(w[keep]^2 * v1[keep]))
  }

  levs <- sort(unique(strata))
  uvs  <- vapply(levs, function(lv) uv_one(strata == lv), numeric(2))
  sum(uvs["U", ]) / sqrt(sum(uvs["V", ]))
}

test_that("stratified FH(0,1) matches the per-stratum pure-R reference", {
  skip_if_not_installed("survival")
  ov  <- survival::ovarian
  jj  <- as.integer(ov$rx != 1)
  z_ref <- ref_swlr_z(ov$futime, ov$fustat, jj, ov$resid.ds,
                      scheme = 0, rho = 0, gamma = 1)
  z     <- as.numeric(survdiff_fast(ov$futime, ov$fustat, ov$rx, 1, side = 1,
                                    weight = "fh", rho = 0, gamma = 1,
                                    strata = ov$resid.ds))
  expect_equal(z, z_ref, tolerance = 1e-8)
})

test_that("stratified mwlrt matches the per-stratum pure-R reference", {
  skip_if_not_installed("survival")
  ov  <- survival::ovarian
  jj  <- as.integer(ov$rx != 1)
  z_ref <- ref_swlr_z(ov$futime, ov$fustat, jj, ov$resid.ds,
                      scheme = 1, t_star = 365)
  z     <- as.numeric(survdiff_fast(ov$futime, ov$fustat, ov$rx, 1, side = 1,
                                    weight = "mwlrt", t_star = 365,
                                    strata = ov$resid.ds))
  expect_equal(z, z_ref, tolerance = 1e-8)
})

test_that("stratified Gehan and Tarone-Ware match the pure-R reference", {
  skip_if_not_installed("survival")
  ov  <- survival::ovarian
  jj  <- as.integer(ov$rx != 1)
  z_gehan_ref <- ref_swlr_z(ov$futime, ov$fustat, jj, ov$resid.ds, scheme = 2)
  z_tw_ref    <- ref_swlr_z(ov$futime, ov$fustat, jj, ov$resid.ds, scheme = 3)
  z_gehan <- as.numeric(survdiff_fast(ov$futime, ov$fustat, ov$rx, 1, side = 1,
                                      weight = "gehan", strata = ov$resid.ds))
  z_tw    <- as.numeric(survdiff_fast(ov$futime, ov$fustat, ov$rx, 1, side = 1,
                                      weight = "tarone-ware",
                                      strata = ov$resid.ds))
  expect_equal(z_gehan, z_gehan_ref, tolerance = 1e-8)
  expect_equal(z_tw,    z_tw_ref,    tolerance = 1e-8)
})

test_that("stratified FH(0,0) reduces to the stratified ordinary log-rank", {
  skip_if_not_installed("survival")
  ov  <- survival::ovarian
  z_swlr <- as.numeric(survdiff_fast(ov$futime, ov$fustat, ov$rx, 1, side = 1,
                                     weight = "fh", rho = 0, gamma = 0,
                                     strata = ov$resid.ds))
  # Stratified ordinary log-rank Z from the unweighted stratified path
  fit_lr <- survdiff_fast(ov$futime, ov$fustat, ov$rx, 1, side = 1,
                          strata = ov$resid.ds)
  expect_equal(z_swlr, as.numeric(fit_lr), tolerance = 1e-8)
})

test_that("a single stratum reduces to the unstratified weighted test", {
  skip_if_not_installed("survival")
  ov  <- survival::ovarian
  one <- rep(1L, nrow(ov))
  z_str <- as.numeric(survdiff_fast(ov$futime, ov$fustat, ov$rx, 1, side = 1,
                                    weight = "fh", rho = 0, gamma = 1,
                                    strata = one))
  z_uns <- as.numeric(survdiff_fast(ov$futime, ov$fustat, ov$rx, 1, side = 1,
                                    weight = "fh", rho = 0, gamma = 1))
  expect_equal(z_str, z_uns, tolerance = 1e-10)
})

test_that("stratified weighted presorted = TRUE matches presorted = FALSE", {
  skip_if_not_installed("survival")
  ov  <- survival::ovarian
  f1  <- survdiff_fast(ov$futime, ov$fustat, ov$rx, 1, side = 1,
                       weight = "fh", rho = 0, gamma = 1, strata = ov$resid.ds)
  ord <- order(ov$resid.ds, ov$futime)
  f2  <- survdiff_fast(ov$futime[ord], ov$fustat[ord], ov$rx[ord], 1, side = 1,
                       presorted = TRUE, weight = "fh", rho = 0, gamma = 1,
                       strata = ov$resid.ds[ord])
  expect_equal(as.numeric(f1), as.numeric(f2), tolerance = 1e-10)
})

test_that("stratified weighted records both weight and strata attributes", {
  skip_if_not_installed("survival")
  ov  <- survival::ovarian
  fit <- survdiff_fast(ov$futime, ov$fustat, ov$rx, 1, side = 1,
                       weight = "fh", rho = 0, gamma = 1, strata = ov$resid.ds)
  expect_equal(attr(fit, "weight"), "fh")
  expect_equal(attr(fit, "strata"), length(unique(ov$resid.ds)))
})

test_that("observed and expected counts include a last subject at risk alone", {
  # Hand calculation: control times 1, 2 and treatment times 3, 4, all events.
  # At t = 4 the only subject at risk is in the treatment group, so O1 = 2 and
  # E1 = 1/2 + 2/3 + 1 + 1 = 19/6, as in survival::survdiff.
  tt <- c(1, 2, 3, 4)
  ee <- c(1, 1, 1, 1)
  gg <- c(0, 0, 1, 1)
  fit <- survdiff_fast(tt, ee, gg, control = 0)
  expect_equal(attr(fit, "O1"), 2)
  expect_equal(attr(fit, "O0"), 2)
  expect_equal(attr(fit, "E1"), 19 / 6, tolerance = 1e-12)
  fw <- survdiff_fast(tt, ee, gg, control = 0, weight = "fh")
  expect_equal(attr(fw, "O1"), 2)
  fs <- survdiff_fast(rep(tt, 2), rep(ee, 2), rep(gg, 2), control = 0,
                      strata = rep(1:2, each = 4))
  expect_equal(attr(fs, "O1"), 4)
  expect_equal(attr(fs, "E1"), 2 * 19 / 6, tolerance = 1e-12)
})

test_that("survdiff_fast validates event coding, group levels, and control", {
  tt <- c(1, 2, 3, 4, 5, 6)
  gg <- c(0, 0, 0, 1, 1, 1)
  expect_error(survdiff_fast(tt, c(1, 2, 1, 2, 1, 2), gg, control = 0),
               "coded as 0")
  expect_error(survdiff_fast(tt, rep(1, 6), c(0, 0, 1, 1, 2, 2), control = 0),
               "two distinct")
  expect_error(survdiff_fast(tt, rep(1, 6), gg, control = 5), "control")
  expect_error(survdiff_fast(tt, rep(1, 6), c(0, NA, 0, 1, 1, 1),
                             control = 0), "missing")
})

test_that("print labels the unweighted stratified log-rank test", {
  skip_if_not_installed("survival")
  ov  <- survival::ovarian
  fit <- survdiff_fast(ov$futime, ov$fustat, ov$rx, 1, strata = ov$resid.ds)
  expect_output(print(fit), "Stratified log-rank test")
})

test_that("survdiff_fast: presorted = TRUE checks the order", {
  time  <- c(3, 1, 2, 4, 6, 5)
  event <- c(1, 1, 0, 1, 1, 0)
  group <- c(0, 1, 0, 1, 0, 1)
  expect_error(survdiff_fast(time, event, group, control = 0, presorted = TRUE), "presorted = FALSE")
  expect_error(survdiff_fast(time, event, group, control = 0, presorted = TRUE, strata = c(1, 1, 1, 2, 2, 2)), "presorted = FALSE")
  expect_error(survdiff_fast(time, event, group, control = 0, presorted = TRUE, weight = "fh", gamma = 1), "presorted = FALSE")
  expect_error(survdiff_fast(time, event, group, control = 0, presorted = TRUE, weight = "fh", gamma = 1, strata = c(1, 1, 1, 2, 2, 2)), "presorted = FALSE")
})

test_that("survdiff_fast: correctly presorted stratified input is accepted", {
  time  <- c(1, 2, 3, 1.5, 2.5, 3.5)
  event <- c(1, 0, 1, 1, 1, 0)
  group <- c(0, 1, 1, 0, 1, 0)
  s     <- c(1, 1, 1, 2, 2, 2)
  expect_equal(
    unclass(survdiff_fast(time, event, group, control = 0, strata = s,
                          presorted = TRUE)),
    unclass(survdiff_fast(time, event, group, control = 0, strata = s)))
  expect_error(survdiff_fast(time, event, group, control = 0, presorted = TRUE,
                             strata = c(1, 2, 1, 2, 1, 2)), "contiguous")
})

test_that("survdiff_fast: mwlrt is not capped when the pooled curve reaches 0 before t_star", {
  # All six subjects have events before t_star = 10, so S(t_star-) = 0 and the
  # weights 1 / S(t-) = 1, 1.2, 1.5, 2, 3, 6 are not capped (as in nphRCT). The
  # expected values were computed independently in Python; the ordinary
  # log-rank Z is -0.6963306, and the stratified Z with a cap of 1 in the first
  # stratum would be 0.5465244.
  tt <- 1:6
  ee <- rep(1, 6)
  gg <- c(0, 1, 0, 1, 0, 1)
  z_mw <- as.numeric(survdiff_fast(tt, ee, gg, control = 0, side = 1,
                                   weight = "mwlrt", t_star = 10))
  expect_equal(z_mw, -0.7734668523, tolerance = 1e-8)
  t2 <- c(2, 4, 6, 8, 10, 12)
  e2 <- c(1, 1, 1, 0, 1, 0)
  g2 <- c(1, 0, 1, 0, 1, 0)
  z_st <- as.numeric(survdiff_fast(c(tt, t2), c(ee, e2), c(gg, g2),
                                   control = 0, side = 1, weight = "mwlrt",
                                   t_star = 10, strata = rep(1:2, each = 6)))
  expect_equal(z_st, 0.0655990630, tolerance = 1e-8)
  # Report whether the external comparison below is run or skipped.
  message("nphRCT available for the mwlrt comparison: ",
          requireNamespace("nphRCT", quietly = TRUE))
  skip_if_not_installed("nphRCT")
  skip_if_not_installed("survival")
  df <- data.frame(time = tt, event = ee, arm = factor(gg))
  ref <- nphRCT::wlrt(survival::Surv(time, event) ~ arm, data = df,
                      method = "mw", t_star = 10)
  expect_equal(abs(z_mw), abs(ref$z), tolerance = 1e-8)
})

test_that("survdiff_fast: a factor event indicator is rejected", {
  expect_error(survdiff_fast(1:6, factor(c(1, 0, 1, 1, 0, 1)), rep(0:1, 3),
                             control = 0), "numeric or logical")
})

test_that("survdiff_fast: the formula interface gives the vector result", {
  set.seed(62)
  n <- 120
  d <- data.frame(tt = round(rexp(n, 0.1), 1), ee = rbinom(n, 1, 0.8),
                  arm = factor(rep(c("ctl", "exp"), n / 2)),
                  s1 = rep(c("a.b", "a"), each = n / 2),
                  s2 = rep(c("c", "b.c", "c"), n / 3))
  ref <- survdiff_fast(d$tt, d$ee, d$arm, control = "ctl")
  expect_identical(survdiff_fast(Surv(tt, ee) ~ arm, data = d), ref)
  expect_identical(survdiff_fast(Surv(tt, ee) ~ arm, d), ref)
  expect_identical(class(ref), c("survdiff_fast", "fastsurv"))
  expect_identical(
    survdiff_fast(Surv(tt, ee) ~ arm, data = d, control = "exp", side = 1,
                  weight = "fh", rho = 0, gamma = 1),
    survdiff_fast(d$tt, d$ee, d$arm, control = "exp", side = 1,
                  weight = "fh", rho = 0, gamma = 1))
  # Strata: the stratum labels only order the strata, so the sums can differ
  # in the last bits from those of other labels
  expect_equal(survdiff_fast(Surv(tt, ee) ~ arm + strata(s1), data = d),
               survdiff_fast(d$tt, d$ee, d$arm, "ctl", strata = d$s1),
               tolerance = 1e-12)
  st12 <- interaction(d$s1, d$s2, sep = "|")
  res12 <- survdiff_fast(Surv(tt, ee) ~ arm + strata(s1, s2), data = d)
  expect_equal(res12, survdiff_fast(d$tt, d$ee, d$arm, "ctl", strata = st12),
               tolerance = 1e-12)
  # Four strata, although pasting the labels with "." would give three
  expect_identical(attr(res12, "strata"), 4L)
  expect_error(survdiff_fast(Surv(tt, ee) ~ arm, data = d, strata = d$s1),
               "strata\\(\\) in the formula")
  expect_error(survdiff_fast(d$tt, d$ee, d$arm, "ctl", data = d),
               "only with a formula")
})

test_that("survdiff_fast: the formula interface passes every argument", {
  set.seed(60)
  n <- 120
  d <- data.frame(tt = round(rexp(n, 0.1), 1) + 0.1, ee = rbinom(n, 1, 0.8),
                  g = rep(c("b", "a"), n / 2), s = rep(c("x", "y", "z"), n / 3))
  d1 <- d[order(d$tt), ]
  d2 <- d[order(d$s, d$tt), ]
  expect_identical(
    survdiff_fast(Surv(tt, ee) ~ g + strata(s), d2, control = "b", side = 1,
                  presorted = TRUE, weight = "mwlrt", t_star = 5),
    survdiff_fast(d2$tt, d2$ee, d2$g, "b", side = 1, presorted = TRUE,
                  strata = d2$s, weight = "mwlrt", t_star = 5))
  expect_identical(
    survdiff_fast(Surv(tt, ee) ~ g, d, control = "b", weight = "fh",
                  rho = 1, gamma = 0.5),
    survdiff_fast(d$tt, d$ee, d$g, "b", weight = "fh", rho = 1, gamma = 0.5))
})
