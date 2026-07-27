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
