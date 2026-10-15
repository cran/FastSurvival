test_that("analysis_fast: required inputs and statistic selection are validated", {
  set.seed(88)
  dat <- simdata_fast(nsim = 5, n = c(40, 40), a.time = c(0, 6), a.prop = 1,
                      e.hazard = list(log(2) / 10, log(2) / 12),
                      d.median = list(30, 30), seed = 88)

  expect_error(analysis_fast(dat, control = 1), "exactly one of")
  expect_error(analysis_fast(dat, control = 1, event.looks = 50,
                             time.looks = 10), "exactly one of")
  expect_error(analysis_fast(dat, control = 1, event.looks = 50,
                             stat = "rmst"), "tau")
  expect_error(analysis_fast(dat, control = 1, event.looks = 50,
                             stat = "km"), "t.eval")
  expect_error(analysis_fast(dat, control = 1, event.looks = 50,
                             stat = "bogus"), "subset of")
  expect_error(analysis_fast(dat, control = 1, event.looks = 50,
                             stat = "logrank", weight = "mwlrt"), "t_star")
})

# Helper: reproduce one (sim, look) cell by applying the administrative cut in
# R and calling the already externally-validated per-statistic wrappers. This
# is the reference the fused kernel must match.
cut_one <- function(dat, sim_id, cutoff) {
  d <- dat[dat$sim == sim_id, ]
  enrolled <- d$accrual_time <= cutoff
  d <- d[enrolled, ]
  before <- (d$accrual_time + d$tte) <= cutoff
  t <- ifelse(before, d$tte, cutoff - d$accrual_time)
  e <- ifelse(before, d$event, 0L)
  ord <- order(t)
  list(time = t[ord], event = as.integer(e[ord]), group = d$group[ord])
}

test_that("analysis_fast logrank/coxph match per-cell wrappers (time looks)", {
  set.seed(101)
  dat <- simdata_fast(nsim = 12, n = c(100, 100), a.time = c(0, 10),
                      a.prop = 1,
                      e.hazard = list(log(2) / 10, log(2) / 14),
                      d.median = list(30, 30), seed = 101)
  looks <- c(12, 22)
  res <- analysis_fast(dat, control = 1, time.looks = looks,
                       stat = c("logrank", "coxph"), side = 2)

  sims <- sort(unique(dat$sim))
  row <- 0L
  for (s in sims) {
    for (cv in looks) {
      row <- row + 1L
      cc <- cut_one(dat, s, cv)
      if (length(cc$time) == 0L) next
      n_ev <- sum(cc$event)
      both <- any(cc$group == 1) && any(cc$group != 1)
      if (n_ev > 0 && both) {
        z_ref <- as.numeric(survdiff_fast(cc$time, cc$event, cc$group,
                                          control = 1, side = 1,
                                          presorted = TRUE))
        expect_equal(res$logrank.z[row], z_ref, tolerance = 1e-10)

        cx <- coxph_fast(cc$time, cc$event, cc$group, control = 1,
                         presorted = TRUE)
        expect_equal(res$cox.coef[row], unname(cx[1]), tolerance = 1e-10)
        expect_equal(res$cox.hr[row],   unname(cx[2]), tolerance = 1e-10)
        expect_equal(res$cox.se[row],   unname(cx[3]), tolerance = 1e-10)
      }
    }
  }
})

test_that("analysis_fast rmst/km/ahsw match per-cell wrappers (event looks)", {
  set.seed(202)
  dat <- simdata_fast(nsim = 12, n = c(120, 120), a.time = c(0, 12),
                      a.prop = 1,
                      e.hazard = list(log(2) / 12, c(log(2) / 12, log(2) / 18)),
                      e.time = list(NULL, c(0, 6, Inf)),
                      d.median = list(36, 36), seed = 202)
  looks <- c(100, 160)
  res <- analysis_fast(dat, control = 1, event.looks = looks,
                       stat = c("rmst", "km", "ahsw"),
                       tau = 12, t.eval = 12, side = 2)

  sims <- sort(unique(dat$sim))
  row <- 0L
  for (s in sims) {
    cal_ev <- with(dat[dat$sim == s & dat$event == 1, ], accrual_time + tte)
    for (cv in looks) {
      row <- row + 1L
      if (cv > length(cal_ev)) next
      cutoff <- sort(cal_ev)[cv]
      cc <- cut_one(dat, s, cutoff)
      both <- any(cc$group == 1) && any(cc$group != 1)
      if (!both) next

      rm <- rmst_fast(cc$time, cc$event, group = cc$group, control = 1,
                      tau = 12, presorted = TRUE)
      expect_equal(res$rmst.diff[row], unname(rm["diff"]), tolerance = 1e-10)
      expect_equal(res$rmst.z[row],    unname(rm["z.diff"]), tolerance = 1e-10)

      is_c <- cc$group == 1
      kc <- survfit_fast(cc$time[is_c], cc$event[is_c], t_eval = 12,
                         presorted = TRUE)
      expect_equal(res$km.surv.ctrl[row], unname(kc["surv"]), tolerance = 1e-10)

      if (sum(cc$event) > 0) {
        ah <- ahsw_fast(cc$time, cc$event, group = cc$group, control = 1,
                        tau = 12, presorted = TRUE)
        expect_equal(res$ahsw.rah[row], unname(ah["rah"]), tolerance = 1e-10)
        expect_equal(res$ahsw.dah[row], unname(ah["dah"]), tolerance = 1e-10)
      }
    }
  }
})

test_that("analysis_fast weighted log-rank matches survdiff_fast (FH, mwlrt)", {
  set.seed(303)
  dat <- simdata_fast(nsim = 12, n = c(110, 110), a.time = c(0, 12),
                      a.prop = 1,
                      e.hazard = list(log(2) / 12, c(log(2) / 12, log(2) / 20)),
                      e.time = list(NULL, c(0, 6, Inf)),
                      d.median = list(36, 36), seed = 303)
  looks <- 150

  for (cfg in list(list(weight = "fh", rho = 0, gamma = 1, t_star = NULL),
                   list(weight = "mwlrt", rho = 0, gamma = 0, t_star = 6))) {
    res <- analysis_fast(dat, control = 1, event.looks = looks,
                         stat = "logrank", weight = cfg$weight,
                         rho = cfg$rho, gamma = cfg$gamma, t_star = cfg$t_star,
                         side = 2)
    sims <- sort(unique(dat$sim))
    row <- 0L
    for (s in sims) {
      row <- row + 1L
      cal_ev <- with(dat[dat$sim == s & dat$event == 1, ], accrual_time + tte)
      if (looks > length(cal_ev)) next
      cutoff <- sort(cal_ev)[looks]
      cc <- cut_one(dat, s, cutoff)
      both <- any(cc$group == 1) && any(cc$group != 1)
      if (sum(cc$event) == 0 || !both) next
      z_ref <- as.numeric(survdiff_fast(
        cc$time, cc$event, cc$group, control = 1, side = 1, presorted = TRUE,
        weight = cfg$weight, rho = cfg$rho, gamma = cfg$gamma,
        t_star = cfg$t_star))
      expect_equal(res$logrank.z[row], z_ref, tolerance = 1e-10,
                   info = cfg$weight)
    }
  }
})

test_that("analysis_fast max-combo matches maxcombo_fast", {
  set.seed(404)
  dat <- simdata_fast(nsim = 10, n = c(120, 120), a.time = c(0, 12),
                      a.prop = 1,
                      e.hazard = list(log(2) / 12, c(log(2) / 12, log(2) / 18)),
                      e.time = list(NULL, c(0, 6, Inf)),
                      d.median = list(36, 36), seed = 404)
  looks <- 150
  res <- analysis_fast(dat, control = 1, event.looks = looks,
                       stat = "maxcombo", side = 1)

  sims <- sort(unique(dat$sim))
  row <- 0L
  for (s in sims) {
    row <- row + 1L
    cal_ev <- with(dat[dat$sim == s & dat$event == 1, ], accrual_time + tte)
    if (looks > length(cal_ev)) next
    cutoff <- sort(cal_ev)[looks]
    cc <- cut_one(dat, s, cutoff)
    both <- any(cc$group == 1) && any(cc$group != 1)
    if (sum(cc$event) == 0 || !both) next
    mc <- maxcombo_fast(cc$time, cc$event, cc$group, control = 1, side = 1,
                        presorted = TRUE)
    expect_equal(res$maxcombo.stat[row], unname(mc["statistic"]),
                 tolerance = 1e-8)
    # The max-combo statistic is deterministic, so it must match tightly. The
    # p-value comes from mvtnorm::pmvnorm (GenzBretz Monte-Carlo integration),
    # which is not reproducible to machine precision, so it is checked loosely.
    expect_equal(res$maxcombo.p[row], unname(mc["p.value"]),
                 tolerance = 1e-2)
  }
})

test_that("analysis_fast by.subgroup produces correct long-form populations", {
  set.seed(505)
  dat <- simdata_fast(nsim = 10, n = c(150, 150), a.time = c(0, 12),
                      a.prop = 1,
                      e.hazard = list(log(2) / 12, log(2) / 16),
                      d.median = list(36, 36),
                      prevalence = c(0.6, 0.4), seed = 505)
  res <- analysis_fast(dat, control = 1, event.looks = 120,
                       stat = "logrank", by.subgroup = TRUE, side = 2)

  expect_true("population" %in% names(res))
  expect_setequal(unique(res$population),
                  c("overall", "subgroup_1", "subgroup_2"))

  # Subgroup rows at a given look share the same cutoff as the overall row.
  ov  <- res[res$population == "overall", ]
  sg1 <- res[res$population == "subgroup_1", ]
  expect_equal(sg1$cutoff, ov$cutoff, tolerance = 1e-10)

  # The subgroup logrank.z matches survdiff_fast on the subgroup subset.
  s1 <- sort(unique(dat$sim))[1]
  cal_ev <- with(dat[dat$sim == s1 & dat$event == 1, ], accrual_time + tte)
  if (120 <= length(cal_ev)) {
    cutoff <- sort(cal_ev)[120]
    d <- dat[dat$sim == s1, ]
    enr <- d$accrual_time <= cutoff
    d <- d[enr, ]
    before <- (d$accrual_time + d$tte) <= cutoff
    t <- ifelse(before, d$tte, cutoff - d$accrual_time)
    e <- ifelse(before, d$event, 0L)
    sel <- d$subgroup == 1
    if (any(sel) && sum(e[sel]) > 0 &&
        any(d$group[sel] == 1) && any(d$group[sel] != 1)) {
      ord <- order(t[sel])
      z_ref <- as.numeric(survdiff_fast(t[sel][ord], as.integer(e[sel])[ord],
                                        d$group[sel][ord], control = 1,
                                        side = 1, presorted = TRUE))
      z_new <- res$logrank.z[res$population == "subgroup_1" &
                               res$sim == s1][1]
      expect_equal(z_new, z_ref, tolerance = 1e-10)
    }
  }
})

test_that("analysis_fast stratified log-rank matches survdiff_fast with strata", {
  set.seed(606)
  dat <- simdata_fast(nsim = 10, n = c(150, 150), a.time = c(0, 12),
                      a.prop = 1,
                      e.hazard = list(log(2) / 12, log(2) / 16),
                      d.median = list(36, 36),
                      prevalence = c(0.5, 0.5), seed = 606)
  looks <- 150
  res <- analysis_fast(dat, control = 1, event.looks = looks,
                       stat = "logrank", strata = "subgroup", side = 2)

  s1 <- sort(unique(dat$sim))[1]
  cal_ev <- with(dat[dat$sim == s1 & dat$event == 1, ], accrual_time + tte)
  if (looks <= length(cal_ev)) {
    cutoff <- sort(cal_ev)[looks]
    d <- dat[dat$sim == s1, ]
    enr <- d$accrual_time <= cutoff
    d <- d[enr, ]
    before <- (d$accrual_time + d$tte) <= cutoff
    t <- ifelse(before, d$tte, cutoff - d$accrual_time)
    e <- as.integer(ifelse(before, d$event, 0L))
    if (sum(e) > 0) {
      z_ref <- as.numeric(survdiff_fast(t, e, d$group,
                                        control = 1, side = 1,
                                        presorted = FALSE,
                                        strata = d$subgroup))
      z_new <- res$logrank.z[res$sim == s1][1]
      expect_equal(z_new, z_ref, tolerance = 1e-10)
    }
  }
})

test_that("analysis_fast medsurv/wkm/wmst match per-cell wrappers (event looks)", {
  set.seed(707)
  dat <- simdata_fast(nsim = 12, n = c(120, 120), a.time = c(0, 12),
                      a.prop = 1,
                      e.hazard = list(log(2) / 12, c(log(2) / 12, log(2) / 18)),
                      e.time = list(NULL, c(0, 6, Inf)),
                      d.median = list(36, 36), seed = 707)
  looks <- c(110, 170)
  res <- analysis_fast(dat, control = 1, event.looks = looks,
                       stat = c("medsurv", "wkm", "wmst"),
                       medsurv.method = "km", wkm.weight = "PF",
                       wmst.tau1 = 2, wmst.tau2 = 12, side = 2)

  sims <- sort(unique(dat$sim))
  row <- 0L
  for (s in sims) {
    cal_ev <- with(dat[dat$sim == s & dat$event == 1, ], accrual_time + tte)
    for (cv in looks) {
      row <- row + 1L
      if (cv > length(cal_ev)) next
      cutoff <- sort(cal_ev)[cv]
      cc <- cut_one(dat, s, cutoff)
      both <- any(cc$group == 1) && any(cc$group != 1)
      if (!both) next
      n_ev <- sum(cc$event)

      # Median survival difference (km variance, default bandwidth).
      md <- medsurv_fast(cc$time, cc$event, group = cc$group, control = 1,
                         method = "km")
      expect_equal(res$medsurv.diff[row], unname(md["diff"]), tolerance = 1e-10)
      expect_equal(res$medsurv.z[row],    unname(md["z"]),    tolerance = 1e-10)

      # Window mean survival time difference over (2, 12].
      wm <- wmst_fast(cc$time, cc$event, group = cc$group, control = 1,
                      tau1 = 2, tau2 = 12)
      expect_equal(res$wmst.diff[row], unname(wm["diff"]), tolerance = 1e-10)
      expect_equal(res$wmst.z[row],    unname(wm["z"]),    tolerance = 1e-10)

      # Weighted Kaplan-Meier (Pepe-Fleming) requires at least one event.
      if (n_ev > 0) {
        wk <- wkm_fast(cc$time, cc$event, cc$group, control = 1, weight = "PF")
        expect_equal(res$wkm.wdiff[row], unname(wk["wdiff"]), tolerance = 1e-10)
        expect_equal(res$wkm.z[row],     unname(wk["z"]),     tolerance = 1e-10)
      }
    }
  }
})

test_that("analysis_fast: new statistic arguments are validated", {
  set.seed(89)
  dat <- simdata_fast(nsim = 5, n = c(40, 40), a.time = c(0, 6), a.prop = 1,
                      e.hazard = list(log(2) / 10, log(2) / 12),
                      d.median = list(30, 30), seed = 89)

  # 'wmst' needs a window through either wmst.tau2 or tau.
  expect_error(analysis_fast(dat, control = 1, event.looks = 50,
                             stat = "wmst"), "wmst.tau2")
  # wmst.tau1 must be non-negative.
  expect_error(analysis_fast(dat, control = 1, event.looks = 50,
                             stat = "wmst", wmst.tau1 = -1, wmst.tau2 = 12),
               "wmst.tau1")
  # medsurv.bw, when supplied, must be positive.
  expect_error(analysis_fast(dat, control = 1, event.looks = 50,
                             stat = "medsurv", medsurv.bw = 0), "medsurv.bw")
  # match.arg guards the method and weight selectors.
  expect_error(analysis_fast(dat, control = 1, event.looks = 50,
                             stat = "medsurv", medsurv.method = "bogus"),
               "should be one of")
  expect_error(analysis_fast(dat, control = 1, event.looks = 50,
                             stat = "wkm", wkm.weight = "bogus"),
               "should be one of")
})

test_that("analysis_fast milestone and ahsw p-values follow side as the wrappers do", {
  dat <- simdata_fast(nsim = 3, n = c(150, 150), a.time = c(0, 12),
                      a.prop = 1,
                      e.hazard = list(log(2) / 10, log(2) / 16), seed = 404)
  # A very late calendar look leaves the data uncut, so each row can be
  # compared with the stand-alone wrappers on the full simulated trial.
  for (mm in c("loglog", "wald", "mover")) {
    res <- analysis_fast(dat, control = 1, time.looks = 1e6,
                         stat = c("milestone", "ahsw"), tau = 12,
                         ms.method = mm, side = 1)
    for (s in 1:3) {
      d  <- dat[dat$sim == s, ]
      ms <- milestone_fast(d$tte, d$event, d$group, control = 1, side = 1,
                           tau = 12, method = mm)
      expect_equal(res$milestone.z[s], ms$statistic, tolerance = 1e-10)
      expect_equal(res$milestone.p[s], ms$p.value, tolerance = 1e-10)
      ah <- ahsw_fast(d$tte, d$event, d$group, control = 1, side = 1,
                      tau = 12)
      expect_equal(res$ahsw.p.rah[s], unname(ah["p.rah"]), tolerance = 1e-10)
      expect_equal(res$ahsw.p.dah[s], unname(ah["p.dah"]), tolerance = 1e-10)
    }
  }
})

test_that("analysis_fast validates group, control, event coding, and looks", {
  dat <- simdata_fast(nsim = 3, n = c(40, 40), a.time = c(0, 6), a.prop = 1,
                      e.hazard = list(log(2) / 10, log(2) / 12), seed = 90)
  expect_error(analysis_fast(dat, control = 3, time.looks = 10), "control")
  dat3 <- dat
  dat3$group[1] <- 3L
  expect_error(analysis_fast(dat3, control = 1, time.looks = 10),
               "two distinct")
  dat_na <- dat
  dat_na$tte[2] <- NA
  expect_error(analysis_fast(dat_na, control = 1, time.looks = 10), "missing")
  dat_ev <- dat
  dat_ev$event[3] <- 2L
  expect_error(analysis_fast(dat_ev, control = 1, time.looks = 10), "coded")
  expect_error(analysis_fast(dat, control = 1, event.looks = 10.5), "whole")
})

test_that("analysis_fast handles factor and character subgroup columns", {
  dat <- simdata_fast(nsim = 4, n = c(100, 100), a.time = c(0, 12),
                      a.prop = 1, e.hazard = list(list(0.10, 0.07), 0.05),
                      prevalence = c(0.5, 0.5), seed = 505)
  ref <- analysis_fast(dat, control = 1, time.looks = 24, by.subgroup = TRUE)
  lab <- c("A", "B")
  dat_f <- dat
  dat_f$subgroup <- factor(lab[dat$subgroup])
  res_f <- analysis_fast(dat_f, control = 1, time.looks = 24,
                         by.subgroup = TRUE)
  expect_equal(res_f$population,
               sub("_2$", "_B", sub("_1$", "_A", ref$population)))
  expect_equal(res_f$n.enrolled, ref$n.enrolled)
  expect_equal(res_f$logrank.z, ref$logrank.z)
  dat_c <- dat
  dat_c$subgroup <- lab[dat$subgroup]
  res_c <- analysis_fast(dat_c, control = 1, time.looks = 24,
                         by.subgroup = TRUE)
  expect_equal(res_c$n.enrolled, ref$n.enrolled)
  expect_equal(res_c$logrank.z, ref$logrank.z)
})

test_that("analysis_fast stratified coxph matches coxph_fast with strata", {
  dat <- simdata_fast(nsim = 4, n = c(150, 150), a.time = c(0, 12),
                      a.prop = 1,
                      e.hazard = list(list(0.10, 0.05), list(0.07, 0.035)),
                      prevalence = c(0.5, 0.5), seed = 707)
  # A very late calendar look leaves the data uncut.
  res <- analysis_fast(dat, control = 1, time.looks = 1e6,
                       stat = c("logrank", "coxph"), strata = "subgroup")
  for (s in 1:4) {
    d  <- dat[dat$sim == s, ]
    cx <- coxph_fast(d$tte, d$event, d$group, control = 1,
                     strata = d$subgroup)
    expect_equal(res$cox.coef[s], unname(cx["coef"]), tolerance = 1e-10)
    expect_equal(res$cox.se[s], unname(cx["se(coef)"]), tolerance = 1e-10)
    lr <- survdiff_fast(d$tte, d$event, d$group, control = 1, side = 1,
                        strata = d$subgroup)
    expect_equal(res$logrank.z[s], as.numeric(lr), tolerance = 1e-10)
  }
})

test_that("analysis_fast cutoff.looks matches per-cell wrappers at per-simulation cutoffs", {
  dat <- simdata_fast(nsim = 10, n = c(100, 100), a.time = c(0, 10),
                      a.prop = 1,
                      e.hazard = list(log(2) / 10, log(2) / 14),
                      d.median = list(30, 30), seed = 303)
  # Simulation-specific cutoffs, with one look that is not reached.
  cut <- cbind(seq(11, 20, length.out = 10), seq(25, 16, length.out = 10))
  cut[3, 2] <- NA
  res <- analysis_fast(dat, control = 1, cutoff.looks = cut,
                       stat = c("logrank", "coxph"))
  expect_equal(nrow(res), 20L)
  expect_true(all(is.na(res$look.value)))

  row <- 0L
  for (s in 1:10) {
    for (l in 1:2) {
      row <- row + 1L
      cv  <- cut[s, l]
      if (is.na(cv)) {
        # Not reached: the statistics of the full (uncut) data.
        expect_false(res$reached[row])
        expect_true(is.na(res$cutoff[row]))
        cc <- cut_one(dat, s, Inf)
      } else {
        expect_true(res$reached[row])
        expect_equal(res$cutoff[row], cv)
        cc <- cut_one(dat, s, cv)
      }
      z_ref <- as.numeric(survdiff_fast(cc$time, cc$event, cc$group,
                                        control = 1, side = 1,
                                        presorted = TRUE))
      expect_equal(res$logrank.z[row], z_ref, tolerance = 1e-10)
      cx <- coxph_fast(cc$time, cc$event, cc$group, control = 1,
                       presorted = TRUE)
      expect_equal(res$cox.coef[row], unname(cx[1]), tolerance = 1e-10)
      expect_equal(res$n.event[row], sum(cc$event))
      expect_equal(res$n.enrolled[row], length(cc$time))
    }
  }
})

test_that("analysis_fast cutoff.looks matches rows by name and validates input", {
  dat <- simdata_fast(nsim = 5, n = c(60, 60), a.time = c(0, 6), a.prop = 1,
                      e.hazard = list(log(2) / 10, log(2) / 12), seed = 404)
  cut <- matrix(c(10, 12, 14, 16, 18), ncol = 1,
                dimnames = list(as.character(5:1), NULL))
  res <- analysis_fast(dat, control = 1, cutoff.looks = cut)
  expect_equal(res$cutoff, c(18, 16, 14, 12, 10))
  vec <- analysis_fast(dat, control = 1, cutoff.looks = c(18, 16, 14, 12, 10))
  expect_equal(vec$logrank.z, res$logrank.z)
  tim <- analysis_fast(dat[dat$sim == 2, ], control = 1, time.looks = 16)
  expect_equal(res$logrank.z[2], tim$logrank.z)

  expect_error(analysis_fast(dat, control = 1, cutoff.looks = matrix(10, 4, 1)),
               "one row per")
  expect_error(analysis_fast(dat, control = 1, cutoff.looks = matrix(-1, 5, 1)),
               "negative")
  expect_error(analysis_fast(dat, control = 1, event.looks = 10,
                             cutoff.looks = matrix(10, 5, 1)),
               "exactly one of")
  bad <- matrix(10, 5, 1, dimnames = list(as.character(2:6), NULL))
  expect_error(analysis_fast(dat, control = 1, cutoff.looks = bad),
               "row names")
})

test_that("analysis_fast max-combo with two or three weights matches maxcombo_fast on both sides", {
  dat <- simdata_fast(nsim = 6, n = c(120, 120), a.time = c(0, 12),
                      a.prop = 1,
                      e.hazard = list(log(2) / 12, c(log(2) / 12, log(2) / 18)),
                      e.time = list(NULL, c(0, 6, Inf)), seed = 405)
  cut <- 1e6
  for (w in list(list(rho = c(0, 0), gamma = c(0, 1)),
                 list(rho = c(0, 0, 1), gamma = c(0, 1, 0)))) {
    for (sd in c(1, 2)) {
      set.seed(1)
      res <- analysis_fast(dat, control = 1, time.looks = cut,
                           stat = "maxcombo", side = sd,
                           mc.rho = w$rho, mc.gamma = w$gamma)
      for (s in 1:6) {
        cc <- cut_one(dat, s, cut)
        mc <- maxcombo_fast(cc$time, cc$event, cc$group, control = 1,
                            side = sd, rho = w$rho, gamma = w$gamma,
                            presorted = TRUE)
        expect_equal(res$maxcombo.stat[s], unname(mc["statistic"]),
                     tolerance = 1e-8)
        # TVPACK (one-sided) is deterministic; the GenzBretz integral
        # (two-sided) carries a randomized integration error, of the
        # order of 1e-4 for singular weight sets (see maxcombo_fast()).
        tol <- if (sd == 1) 1e-8 else 1e-3
        expect_lt(abs(res$maxcombo.p[s] - unname(mc["p.value"])), tol)
      }
    }
  }
})

test_that("analysis_fast max-combo with mc.alpha keeps the decisions and bounds the p-values", {
  dat <- simdata_fast(nsim = 60, n = c(100, 100), a.time = c(0, 12),
                      a.prop = 1,
                      e.hazard = list(log(2) / 12, c(log(2) / 12, log(2) / 20)),
                      e.time = list(NULL, c(0, 4, Inf)), seed = 406)
  lev <- c(0.01, 0.025)
  for (sd in c(1, 2)) {
    for (w in list(list(rho = c(0, 0), gamma = c(0, 1)),
                   list(rho = c(0, 0, 1, 1), gamma = c(0, 1, 0, 1)))) {
      K <- length(w$rho)
      set.seed(2)
      ex <- analysis_fast(dat, control = 1, time.looks = c(15, 30),
                          stat = "maxcombo", side = sd,
                          mc.rho = w$rho, mc.gamma = w$gamma)
      set.seed(2)
      bd <- analysis_fast(dat, control = 1, time.looks = c(15, 30),
                          stat = "maxcombo", side = sd,
                          mc.rho = w$rho, mc.gamma = w$gamma, mc.alpha = lev)
      expect_false("maxcombo.p.exact" %in% names(ex))
      expect_true(is.logical(bd$maxcombo.p.exact))
      expect_equal(bd$maxcombo.stat, ex$maxcombo.stat)
      ok <- !is.na(ex$maxcombo.p)
      expect_true(all(ok))
      a    <- lev[bd$look]
      p_lo <- if (sd == 1) pnorm(bd$maxcombo.stat) else
        2 * pnorm(-bd$maxcombo.stat)
      p_hi <- pmin(1, K * p_lo)
      # The exact p-values lie within the Bonferroni bounds (up to the
      # integration error).
      expect_true(all(ex$maxcombo.p >= p_lo - 1e-4 &
                        ex$maxcombo.p <= p_hi + 1e-4))
      ix <- bd$maxcombo.p.exact
      expect_true(any(!ix))
      # Rows that are integrated agree with the run without mc.alpha.
      tol <- if (sd == 1 && K <= 3) 1e-10 else 1e-3
      expect_lt(max(c(0, abs(bd$maxcombo.p[ix] - ex$maxcombo.p[ix]))), tol)
      # Rows that are not integrated report the bound on the side of the level.
      expect_equal(bd$maxcombo.p[!ix],
                   ifelse(p_lo[!ix] > a[!ix], p_lo[!ix], p_hi[!ix]))
      expect_true(all(ix | p_lo > a | p_hi <= a))
      # The decisions at the levels are those of the exact p-values (rows
      # within the integration error of the level are excluded).
      clear <- abs(ex$maxcombo.p - a) > 1e-3
      expect_equal((bd$maxcombo.p <= a)[clear], (ex$maxcombo.p <= a)[clear])
    }
  }
  expect_error(analysis_fast(dat, control = 1, time.looks = c(15, 30),
                             stat = "maxcombo",
                             mc.alpha = c(0.01, 0.02, 0.03)),
               "mc.alpha")
  expect_error(analysis_fast(dat, control = 1, time.looks = 30,
                             stat = "maxcombo", mc.alpha = 1),
               "mc.alpha")
})

test_that("analysis_fast: an event target beyond the integer range is not reached", {
  df <- simdata_fast(nsim = 3, n = c(30, 30), a.time = c(0, 6), a.rate = 10,
                     e.median = list(12, 18), seed = 1)
  res <- analysis_fast(df, control = 1, event.looks = 3e9)
  expect_true(all(!res$reached))
  expect_true(all(is.na(res$cutoff)))
  expect_equal(res$look.value, rep(3e9, 3))
})

test_that("analysis_fast: wkm is NA at an unreached look with infinite times", {
  # Cure model without dropout: subjects without an event have tte = Inf, so
  # at a look that is not reached their observed time is infinite and the
  # integral of the weighted Kaplan-Meier test has no finite upper limit.
  d <- simdata_fast(nsim = 5, n = c(100, 100), a.time = c(0, 12),
                    a.rate = 200 / 12, e.hazard = list(c(0.1, 0), c(0.1, 0)),
                    e.time = c(0, 24, Inf), seed = 5)
  expect_true(any(is.infinite(d$tte)))
  r <- analysis_fast(d, control = 1, event.looks = c(30, 1000),
                     stat = c("logrank", "wkm"))
  un <- r[!r$reached, ]
  expect_equal(nrow(un), 5L)
  # The weighted difference comes from the C++ core without arithmetic, so
  # it is NA and not NaN; the derived columns are NA.
  expect_true(all(is.na(un[["wkm.wdiff"]]) & !is.nan(un[["wkm.wdiff"]])))
  for (cn in c("wkm.lower", "wkm.upper", "wkm.z", "wkm.p")) {
    expect_true(all(is.na(un[[cn]])), info = cn)
  }
  expect_true(all(is.finite(un[["logrank.z"]])))
  # Reached looks (finite observed times) are unaffected.
  rc <- r[r$reached, ]
  expect_equal(nrow(rc), 5L)
  expect_true(all(is.finite(rc[["wkm.wdiff"]]) & is.finite(rc[["wkm.z"]])))
})

test_that("analysis_fast: medsurv is NA when a median at 0.5 runs up to Inf", {
  # The look is not reached, so the full data are used. In the treatment group
  # the Kaplan-Meier estimate is exactly 0.5 after time 5 and the remaining
  # subjects have tte = Inf, so its median is not defined (see medsurv_fast).
  d <- data.frame(sim = 1L, group = rep(1:2, each = 10), accrual_time = 0,
                  tte = c(1:10, 1:5, rep(Inf, 5)),
                  event = c(rep(1L, 10), rep(1L, 5), rep(0L, 5)))
  for (m in c("km", "nph")) {
    r <- analysis_fast(d, control = 1, event.looks = 1000, stat = "medsurv",
                       medsurv.method = m)
    expect_false(r$reached)
    expect_equal(r$medsurv.ctrl, 5.5)
    for (cn in c("medsurv.trt", "medsurv.diff", "medsurv.diff.lower",
                 "medsurv.diff.upper", "medsurv.z", "medsurv.p")) {
      expect_true(is.na(r[[cn]]), info = paste(m, cn))
    }
  }
})

test_that("analysis_fast: mwlrt is not capped when the pooled curve reaches 0 before t_star", {
  # Same data as in test-survdiff_fast.R; the expected values were computed
  # independently in Python (see survdiff_fast).
  d <- data.frame(sim = 1L, group = c(0, 1, 0, 1, 0, 1), accrual_time = 0,
                  tte = 1:6, event = 1L)
  r <- analysis_fast(d, control = 0, time.looks = 100, weight = "mwlrt",
                     t_star = 10)
  expect_equal(r$logrank.z, -0.7734668523, tolerance = 1e-8)
  d2 <- rbind(cbind(d, s = 1L),
              data.frame(sim = 1L, group = c(1, 0, 1, 0, 1, 0),
                         accrual_time = 0, tte = c(2, 4, 6, 8, 10, 12),
                         event = c(1L, 1L, 1L, 0L, 1L, 0L), s = 2L))
  r2 <- analysis_fast(d2, control = 0, time.looks = 100, weight = "mwlrt",
                      t_star = 10, strata = "s")
  expect_equal(r2$logrank.z, 0.0655990630, tolerance = 1e-8)
})

test_that("analysis_fast: tte must be non-negative and -0 is treated as 0", {
  d <- data.frame(sim = 1L, group = rep(0:1, each = 4), accrual_time = 0,
                  tte = c(-0, 1, 2, 4, 0.5, 1.5, 3, 5), event = 1L)
  # The first time is a negative zero.
  expect_identical(1 / d$tte[1], -Inf)
  r <- analysis_fast(d, control = 0, time.looks = 100)
  z_ref <- as.numeric(survdiff_fast(d$tte, d$event, d$group, control = 0,
                                    side = 1))
  expect_equal(r$logrank.z, z_ref, tolerance = 1e-12)
  # Python, with the event at time 0 first: -0.6414478 (+0.6414478 when the
  # negative zero is sorted last).
  expect_equal(r$logrank.z, -0.6414478072, tolerance = 1e-8)
  d_neg <- d
  d_neg$tte[1] <- -0.3
  expect_error(analysis_fast(d_neg, control = 0, time.looks = 100),
               "non-negative")
  d_chr <- d
  d_chr$tte <- as.character(d_chr$tte)
  expect_error(analysis_fast(d_chr, control = 0, time.looks = 100),
               "numeric and non-negative")
  d_fac <- d
  d_fac$event <- factor(d_fac$event)
  expect_error(analysis_fast(d_fac, control = 0, time.looks = 100),
               "numeric or logical")
})

test_that("analysis_fast: no warning when tau exceeds the follow-up at a look", {
  d <- simdata_fast(nsim = 1, n = c(100, 100), a.time = c(0, 12), a.prop = 1,
                    e.median = list(12, 18), seed = 3)
  expect_warning(
    r <- analysis_fast(d, control = 1, time.looks = 10, stat = "rmst",
                       tau = 24),
    regexp = NA)
  inc <- d$accrual_time <= 10
  obs <- pmin(d$tte, 10 - d$accrual_time)[inc]
  evt <- as.integer(d$event == 1 & d$accrual_time + d$tte <= 10)[inc]
  # The stand-alone function warns, and gives the same values.
  expect_warning(f <- rmst_fast(obs, evt, d$group[inc], control = 1, tau = 24),
                 "carried forward")
  expect_equal(r$rmst.ctrl, unname(f["rmst.ctrl"]), tolerance = 1e-10)
  expect_equal(r$rmst.trt, unname(f["rmst.trt"]), tolerance = 1e-10)
})

test_that("analysis_fast: conf.level applies to every confidence interval", {
  dat <- simdata_fast(nsim = 2, n = c(150, 150), a.time = c(0, 12),
                      a.prop = 1,
                      e.hazard = list(log(2) / 10, log(2) / 16), seed = 404)
  res <- analysis_fast(dat, control = 1, time.looks = 1e6,
                       stat = c("milestone", "wmst"), tau = 12,
                       conf.level = 0.9)
  for (s in 1:2) {
    d  <- dat[dat$sim == s, ]
    ms <- milestone_fast(d$tte, d$event, d$group, control = 1, tau = 12,
                         conf.level = 0.9)
    expect_equal(res$milestone.diff.lower[s], ms$diff.lower, tolerance = 1e-10)
    expect_equal(res$milestone.diff.upper[s], ms$diff.upper, tolerance = 1e-10)
    wm <- wmst_fast(d$tte, d$event, d$group, control = 1, tau2 = 12,
                    conf.level = 0.9)
    expect_equal(res$wmst.diff.lower[s], unname(wm["lower.diff"]),
                 tolerance = 1e-10)
    expect_equal(res$wmst.diff.upper[s], unname(wm["upper.diff"]),
                 tolerance = 1e-10)
  }
})

test_that("analysis_fast: ahsw reports finite average hazards as ahsw_fast does", {
  # No treatment-group event up to tau = 5, so the treatment average hazard is
  # 0 and the contrasts are not defined.
  d <- data.frame(sim = 1L, group = rep(1:2, each = 6), accrual_time = 0,
                  tte = c(1:6, 6.5:11.5), event = 1L)
  ah <- ahsw_fast(d$tte, d$event, d$group, control = 1, tau = 5)
  r <- analysis_fast(d, control = 1, time.looks = 100, stat = "ahsw", tau = 5)
  expect_equal(unname(ah["ah.trt"]), 0)
  expect_equal(r$ahsw.ah.ctrl, unname(ah["ah.ctrl"]), tolerance = 1e-12)
  expect_equal(r$ahsw.ah.trt, 0)
  expect_true(is.na(r$ahsw.rah) && is.na(r$ahsw.dah))
})

test_that("analysis_fast: several strata columns do not merge on their labels", {
  # The pasted labels of ("a.b", "c") and ("a", "b.c") are both "a.b.c".
  # Python (independent): Z = -0.8839063432 with the three strata, and
  # -0.5590990412 if the first two were merged.
  d <- data.frame(
    sim = 1L, accrual_time = 0,
    group = rep(rep(0:1, each = 4), 3),
    tte = c(1, 3, 5, 7, 2, 4, 6, 8, 10, 12, 14, 16, 11, 13, 15, 17,
            2.5, 4.5, 6.5, 8.5, 3.5, 5.5, 7.5, 9.5),
    event = c(1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0,
              1, 0, 1, 1, 1, 1, 0, 1),
    s1 = rep(c("a.b", "a", "a"), each = 8),
    s2 = rep(c("c", "b.c", "c"), each = 8),
    stringsAsFactors = FALSE)
  r <- analysis_fast(d, control = 0, time.looks = 100,
                     strata = c("s1", "s2"))
  z_ref <- as.numeric(survdiff_fast(d$tte, d$event, d$group, control = 0,
                                    side = 1,
                                    strata = paste(d$s1, d$s2, sep = "|")))
  expect_equal(r$logrank.z, z_ref, tolerance = 1e-12)
  expect_equal(r$logrank.z, -0.8839063432, tolerance = 1e-8)
})

test_that("analysis_fast: further arguments are validated", {
  dat <- simdata_fast(nsim = 3, n = c(40, 40), a.time = c(0, 6), a.prop = 1,
                      e.hazard = list(log(2) / 10, log(2) / 12), seed = 91)
  expect_error(analysis_fast(dat, control = 1, time.looks = 10,
                             stat = "maxcombo", mc.rho = numeric(0),
                             mc.gamma = numeric(0)),
               "at least one Fleming-Harrington weight")
  expect_error(analysis_fast(dat, control = 1, time.looks = -1),
               "'time.looks' must be positive")
  expect_error(analysis_fast(dat, control = 1, event.looks = Inf),
               "'event.looks' must be positive")
  expect_error(analysis_fast(dat, control = 1, time.looks = 10,
                             conf.level = NA), "conf.level")
  expect_error(analysis_fast(dat, control = 1, time.looks = 10,
                             conf.level = c(0.9, 0.95)), "conf.level")
  expect_error(analysis_fast(dat, control = 1, time.looks = 10,
                             by.subgroup = "yes"), "by.subgroup")
})
