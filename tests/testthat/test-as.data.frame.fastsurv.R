# Two groups with unsorted and tied times, shared by the tests of this file
fs_data <- function(seed) {
  set.seed(seed)
  n <- 150
  data.frame(tt = round(rexp(n, 0.1), 1), ee = rbinom(n, 1, 0.8),
             g = rep(c(0, 1), n / 2))
}

test_that("as.data.frame.fastsurv gives the same columns for every result", {
  d   <- fs_data(81)
  ord <- order(d$tt)
  fits <- list(
    survfit   = survfit_fast(d$tt[ord], d$ee[ord], t_eval = 5),
    survdiff  = survdiff_fast(d$tt, d$ee, d$g, control = 0),
    coxph     = coxph_fast(d$tt, d$ee, d$g, control = 0),
    rmst1     = rmst_fast(d$tt, d$ee, tau = 5),
    rmst2     = rmst_fast(d$tt, d$ee, d$g, control = 0, tau = 5),
    wmst1     = wmst_fast(d$tt, d$ee, tau1 = 1, tau2 = 5),
    wmst2     = wmst_fast(d$tt, d$ee, d$g, control = 0, tau1 = 1, tau2 = 5),
    med1      = medsurv_fast(d$tt, d$ee),
    med2      = medsurv_fast(d$tt, d$ee, d$g, control = 0),
    milestone = milestone_fast(d$tt, d$ee, d$g, control = 0, tau = 5),
    maxcombo  = maxcombo_fast(d$tt, d$ee, d$g, control = 0,
                              rho = c(0, 0), gamma = c(0, 1)),
    rmw       = rmw_fast(d$tt, d$ee, d$g, control = 0),
    wkm       = wkm_fast(d$tt, d$ee, d$g, control = 0),
    ahsw      = ahsw_fast(d$tt, d$ee, d$g, control = 0, tau = 5),
    ahr       = ahr_fast(d$tt, d$ee, d$g, control = 0, tau = 5)
  )
  terms <- list(
    survfit = "surv", survdiff = "logrank", coxph = "coef", rmst1 = "rmst",
    rmst2 = c("rmst.ctrl", "rmst.trt", "diff", "log.ratio"),
    wmst1 = "wmst", wmst2 = c("wmst.control", "wmst.treatment", "diff"),
    med1 = "median", med2 = c("median.control", "median.treatment", "diff"),
    milestone = c("surv.control", "surv.treatment", "diff"),
    maxcombo = c("maxcombo", "FH(0,0)", "FH(0,1)"),
    rmw = c("rmw", "logrank", "mwlrt"), wkm = "wdiff",
    ahsw = c("ah.ctrl", "ah.trt", "log.rah", "dah"),
    ahr = c("log.ahr", "theta")
  )
  cols <- c("term", "estimate", "std.error", "conf.low", "conf.high",
            "statistic", "p.value")
  for (nm in names(fits)) {
    expect_true(inherits(fits[[nm]], "fastsurv"), info = nm)
    tab <- as.data.frame(fits[[nm]])
    expect_identical(names(tab), cols, info = nm)
    expect_identical(tab$term, terms[[nm]], info = nm)
    expect_null(attr(tab, "estimates"))
    expect_null(attr(tab, "conf.level"))
  }
  # Results stack with rbind(), and data.frame() gives the same table
  st <- do.call(rbind, lapply(fits, as.data.frame))
  expect_identical(nrow(st), sum(lengths(terms)))
  expect_equal(data.frame(fits$coxph), as.data.frame(fits$coxph))
  expect_identical(
    rownames(as.data.frame(fits$rmst2, row.names = letters[1:4])),
    letters[1:4])
})

test_that("as.data.frame.fastsurv reports the elements of the estimates", {
  d   <- fs_data(82)
  ord <- order(d$tt)

  sf <- survfit_fast(d$tt[ord], d$ee[ord], t_eval = 5)
  t1 <- as.data.frame(sf)
  expect_equal(unlist(t1[1, c("estimate", "std.error", "conf.low",
                              "conf.high")]),
               unclass(sf)[c("surv", "std.err", "lower", "upper")],
               ignore_attr = TRUE)
  expect_true(is.na(t1$statistic) && is.na(t1$p.value))

  # Log hazard ratio: limits coef -/+ q se at conf.level = 0.9 and the
  # one-sided p-value in the lower tail
  cx  <- coxph_fast(d$tt, d$ee, d$g, control = 0, side = 1, conf.level = 0.9)
  tcx <- as.data.frame(cx)
  b   <- cx[["coef"]]
  se  <- cx[["se(coef)"]]
  expect_equal(tcx$estimate, b)
  expect_equal(tcx$std.error, se)
  expect_equal(tcx$conf.low, b - qnorm(0.95) * se, tolerance = 1e-12)
  expect_equal(tcx$conf.high, b + qnorm(0.95) * se, tolerance = 1e-12)
  expect_equal(tcx$statistic, b / se)
  expect_equal(tcx$p.value, pnorm(b / se))
  expect_equal(as.data.frame(coxph_fast(d$tt, d$ee, d$g, 0))$p.value,
               2 * pnorm(-abs(b / se)))

  # Window mean and median survival: per-group rows and the difference
  for (fit in list(wmst_fast(d$tt, d$ee, d$g, 0, tau1 = 1, tau2 = 5),
                   medsurv_fast(d$tt, d$ee, d$g, 0, side = 1))) {
    tb <- as.data.frame(fit)
    v  <- unclass(fit)
    stem <- if (inherits(fit, "wmst_fast")) "wmst" else "median"
    expect_equal(tb$estimate,
                 unname(v[c(paste0(stem, c(".control", ".treatment")),
                            "diff")]))
    expect_equal(tb$std.error,
                 unname(v[c("se.control", "se.treatment", "se.diff")]))
    expect_equal(tb$conf.low,
                 unname(v[c("lower.control", "lower.treatment",
                            "lower.diff")]))
    expect_equal(tb$conf.high,
                 unname(v[c("upper.control", "upper.treatment",
                            "upper.diff")]))
    expect_equal(tb$statistic, c(NA, NA, unname(v["z"])))
    expect_equal(tb$p.value, c(NA, NA, unname(v["p"])))
  }

  wk <- wkm_fast(d$tt, d$ee, d$g, 0)
  expect_equal(unlist(as.data.frame(wk)[1, -1L]),
               unclass(wk)[c("wdiff", "se", "lower", "upper", "z", "p")],
               ignore_attr = TRUE)
})

test_that("as.data.frame.fastsurv: ratios on the log scale and derived SEs", {
  d <- fs_data(83)
  q <- qnorm(0.975)

  # Ratio of restricted mean survival times: the standard error recovered
  # from the limits reproduces the stored z statistic and p-value
  rm <- rmst_fast(d$tt, d$ee, d$g, control = 0, tau = 5)
  lr <- as.data.frame(rm)[4L, ]
  expect_equal(exp(lr$estimate), rm[["ratio"]])
  expect_equal(exp(c(lr$conf.low, lr$conf.high)),
               c(rm[["ratio.lower"]], rm[["ratio.upper"]]))
  expect_equal(lr$estimate - q * lr$std.error, lr$conf.low, tolerance = 1e-12)
  expect_equal(lr$estimate / lr$std.error, rm[["z.ratio"]], tolerance = 1e-10)
  expect_equal(lr$p.value, 2 * pnorm(-abs(lr$statistic)), tolerance = 1e-10)
  df <- as.data.frame(rm)[3L, ]
  expect_equal(c(df$estimate, df$std.error, df$statistic, df$p.value),
               c(rm[["diff"]], rm[["se.diff"]], rm[["z.diff"]],
                 rm[["p.diff"]]))

  # Average hazard with survival weight, two-sided and one-sided: the
  # derived statistics give the stored p-values
  for (sd in 1:2) {
    ah <- ahsw_fast(d$tt, d$ee, d$g, control = 0, tau = 5, side = sd)
    ta <- as.data.frame(ah)
    expect_equal(ta$estimate[1:2], c(ah[["ah.ctrl"]], ah[["ah.trt"]]))
    expect_equal(exp(ta$estimate[3L]), ah[["rah"]])
    expect_equal(exp(c(ta$conf.low[3L], ta$conf.high[3L])),
                 c(ah[["rah.lower"]], ah[["rah.upper"]]))
    expect_equal(c(ta$conf.low[4L], ta$conf.high[4L]),
                 c(ah[["dah.lower"]], ah[["dah.upper"]]))
    p_from_z <- if (sd == 1L) {
      pnorm(ta$statistic[3:4])
    } else {
      2 * pnorm(-abs(ta$statistic[3:4]))
    }
    expect_equal(p_from_z, c(ah[["p.rah"]], ah[["p.dah"]]),
                 tolerance = 1e-10)
    expect_equal(ta$p.value[3:4], c(ah[["p.rah"]], ah[["p.dah"]]))
  }

  # Kalbfleisch-Prentice average hazard ratio: the log-scale row and the
  # share of the treatment group with the primary test
  ar <- ahr_fast(d$tt, d$ee, d$g, control = 0, tau = 5)
  tr <- as.data.frame(ar)
  expect_equal(exp(tr$estimate[1L]), ar$ahr)
  expect_equal(exp(c(tr$conf.low[1L], tr$conf.high[1L])),
               c(ar$lower, ar$upper))
  expect_equal(c(tr$std.error[1L], tr$statistic[1L], tr$p.value[1L]),
               c(ar$se.loghr, ar$z.loghr, ar$p.value.loghr))
  expect_equal(tr$estimate[2L], unname(ar$theta[2L]))
  # With null.ahr = 1 the primary statistic is (theta - 1/2) / se.theta
  expect_equal(tr$statistic[2L], (tr$estimate[2L] - 0.5) / tr$std.error[2L],
               tolerance = 1e-12)
  expect_equal(tr$p.value[2L], ar$p.value)
  expect_true(all(is.na(c(tr$conf.low[2L], tr$conf.high[2L]))))
})

test_that("as.data.frame.fastsurv: milestone survival by method", {
  d <- fs_data(84)
  q <- qnorm(0.975)
  for (m in c("wald", "loglog", "mover")) {
    ms <- milestone_fast(d$tt, d$ee, d$g, control = 0, tau = 5, method = m)
    tm <- as.data.frame(ms)
    expect_equal(tm$estimate,
                 c(ms$surv[["control"]], ms$surv[["treatment"]], ms$diff))
    expect_equal(tm$std.error[1:2], unname(ms$std.err))
    expect_equal(tm$conf.low, c(unname(ms$surv.lower), ms$diff.lower))
    expect_equal(tm$conf.high, c(unname(ms$surv.upper), ms$diff.upper))
    expect_equal(tm$statistic[3L], ms$statistic)
    expect_equal(tm$p.value[3L], ms$p.value)
    if (m == "wald") {
      # The Wald interval is the difference -/+ q times the standard error
      expect_equal(tm$estimate[3L] + c(-1, 1) * q * tm$std.error[3L],
                   c(ms$diff.lower, ms$diff.upper), tolerance = 1e-12)
    } else {
      expect_true(is.na(tm$std.error[3L]))
    }
  }
})

test_that("as.data.frame.fastsurv: tests carry the statistic and p-value", {
  d <- fs_data(85)
  # One-sided log-rank: the p-value is the lower tail of Z, computed here
  # from the observed and expected counts
  s1 <- survdiff_fast(d$tt, d$ee, d$g, control = 0, side = 1)
  z  <- (attr(s1, "O1") - attr(s1, "E1")) / sqrt(attr(s1, "V1"))
  t1 <- as.data.frame(s1)
  expect_equal(t1$statistic, z, tolerance = 1e-12)
  expect_equal(t1$p.value, pnorm(z), tolerance = 1e-12)
  expect_true(all(is.na(unlist(t1[, c("estimate", "std.error", "conf.low",
                                      "conf.high")]))))
  expect_identical(
    as.data.frame(survdiff_fast(d$tt, d$ee, d$g, 0, weight = "fh",
                                gamma = 1))$term, "fh")

  mc <- maxcombo_fast(d$tt, d$ee, d$g, control = 0, rho = c(0, 0, 1),
                      gamma = c(0, 1, 0), side = 1)
  tm <- as.data.frame(mc)
  expect_equal(tm$statistic, c(mc[["statistic"]], unname(attr(mc, "z"))))
  expect_equal(tm$p.value, c(mc[["p.value"]], NA, NA, NA))
  expect_identical(tm$term, c("maxcombo", names(attr(mc, "z"))))

  rw <- rmw_fast(d$tt, d$ee, d$g, control = 0, side = 1)
  tw <- as.data.frame(rw)
  expect_equal(tw$statistic, c(rw[["statistic"]], attr(rw, "z")[["logrank"]],
                               attr(rw, "z")[["mwlrt"]]))
  expect_equal(tw$p.value, c(rw[["p.value"]], NA, NA))
})

test_that("as.data.frame.fastsurv: two-sided log-rank p-value of survival", {
  d <- fs_data(86)
  t2 <- as.data.frame(survdiff_fast(d$tt, d$ee, d$g, control = 0))
  # Report whether the external comparison below is run or skipped.
  message("survival available for the survdiff p-value comparison: ",
          requireNamespace("survival", quietly = TRUE))
  skip_if_not_installed("survival")
  ref <- survival::survdiff(survival::Surv(tt, ee) ~ g, data = d)
  expect_equal(t2$statistic, ref$chisq, tolerance = 1e-10)
  expect_equal(t2$p.value, pchisq(ref$chisq, df = 1, lower.tail = FALSE),
               tolerance = 1e-10)
})

test_that("as.data.frame.fastsurv keeps a row for an unavailable estimate", {
  na_cx <- coxph_fast(c(1, 2, 3, 4), c(0, 0, 0, 0), c(0, 0, 1, 1),
                      control = 0)
  tab <- as.data.frame(na_cx)
  expect_identical(tab$term, "coef")
  expect_true(all(is.na(unlist(tab[, -1L]))))
})
