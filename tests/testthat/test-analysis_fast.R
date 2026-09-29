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
                             stat = "medsurv", medsurv.method = "bogus"))
  expect_error(analysis_fast(dat, control = 1, event.looks = 50,
                             stat = "wkm", wkm.weight = "bogus"))
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
