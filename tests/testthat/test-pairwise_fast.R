make_karm <- function(nsim = 100, seed = 808) {
  simdata_fast(
    nsim     = nsim,
    n        = c(120, 120, 120),
    a.time   = c(0, 12),
    a.rate   = 360 / 12,
    e.median = list(12, 16, 20),
    seed     = seed
  )
}

test_that("pairwise_fast: fixed time.looks produces one block per contrast", {
  nsim <- 100
  df   <- make_karm(nsim)
  pw <- pairwise_fast(df, control = 1, time.looks = 30,
                      stat = "logrank", side = 1)

  # One row per arm per look per simulation; two experimental arms, one look.
  expect_equal(sort(unique(pw$arm)), c(2, 3))
  expect_equal(nrow(pw), 2 * nsim)
  expect_true(all(c("arm", "sim", "look", "cutoff", "reached",
                    "logrank.z", "logrank.p") %in% names(pw)))
  expect_true(all(is.finite(pw$logrank.z)))
  expect_true(all(pw$logrank.p >= 0 & pw$logrank.p <= 1))
})

test_that("pairwise_fast: fixed-time contrast equals a direct analysis_fast call", {
  df <- make_karm(100)
  pw <- pairwise_fast(df, control = 1, time.looks = 30,
                      stat = "logrank", side = 1)

  sub2 <- df[df$group %in% c(1, 2), ]
  ref  <- analysis_fast(sub2, control = 1, time.looks = 30,
                        stat = "logrank", side = 1)

  a2 <- pw[pw$arm == 2, ]
  a2 <- a2[match(ref$sim, a2$sim), ]
  expect_equal(a2$logrank.z, ref$logrank.z, tolerance = 1e-10)
  expect_equal(a2$logrank.p, ref$logrank.p, tolerance = 1e-10)
})

test_that("pairwise_fast: Bonferroni multiplies p by the number of contrasts", {
  df <- make_karm(100)
  pw <- pairwise_fast(df, control = 1, time.looks = 30,
                      stat = "logrank", side = 1, adjust = "bonferroni")
  expect_true("p.adj" %in% names(pw))
  expect_equal(pw$p.adj, pmin(1, pw$logrank.p * 2), tolerance = 1e-12)
})

test_that("pairwise_fast: 'arms' restricts the set of contrasts", {
  df <- make_karm(60)
  pw <- pairwise_fast(df, control = 1, arms = 2, time.looks = 30,
                      stat = "logrank", side = 1)
  expect_equal(unique(pw$arm), 2)
})

test_that("pairwise_fast: event-driven mode shares the primary cutoff", {
  nsim <- 100
  df   <- make_karm(nsim)
  pw <- pairwise_fast(df, control = 1, event.looks = 200, primary = 3,
                      stat = "logrank", side = 1)

  expect_equal(sort(unique(pw$arm)), c(2, 3))
  expect_equal(nrow(pw), 2 * nsim)
  # The event target is reached in every simulation here (no dropout).
  expect_true(all(pw$reached))
  # Both contrasts at a given simulation use the same calendar cutoff.
  cut2 <- pw$cutoff[pw$arm == 2]
  cut3 <- pw$cutoff[pw$arm == 3]
  expect_equal(cut2, cut3, tolerance = 1e-10)

  # The primary contrast reproduces a direct event-driven analysis.
  sub3 <- df[df$group %in% c(1, 3), ]
  ref3 <- analysis_fast(sub3, control = 1, event.looks = 200,
                        stat = "logrank", side = 1)
  a3 <- pw[pw$arm == 3, ]
  a3 <- a3[match(ref3$sim, a3$sim), ]
  expect_equal(a3$logrank.z, ref3$logrank.z, tolerance = 1e-8)
  expect_equal(a3$cutoff,    ref3$cutoff,    tolerance = 1e-8)
})

test_that("pairwise_fast: input validation", {
  df <- make_karm(30)

  expect_error(pairwise_fast(df, control = 1),
               "exactly one of 'event.looks', 'time.looks', or 'cutoff.looks'")
  expect_error(pairwise_fast(df, control = 1, event.looks = 100, time.looks = 30),
               "exactly one of 'event.looks', 'time.looks', or 'cutoff.looks'")
  expect_error(pairwise_fast(df, control = 1, event.looks = 100),
               "'primary' must name")
  expect_error(pairwise_fast(df, control = 9, time.looks = 30),
               "not among the group labels")
  expect_error(pairwise_fast(df, control = 1, arms = c(1, 2), time.looks = 30),
               "must not include the control")
  expect_error(pairwise_fast(df, control = 1, time.looks = 30, by.subgroup = TRUE),
               "does not support 'by.subgroup'")
})

test_that("pairwise_fast: event-driven dropout and pipeline counts match analysis_fast", {
  df <- simdata_fast(nsim = 30, n = c(120, 120, 120), a.time = c(0, 12),
                     a.rate = 360 / 12, e.median = list(12, 16, 20),
                     d.hazard = 0.02, seed = 909)
  pw <- pairwise_fast(df, control = 1, event.looks = 150, primary = 3,
                      stat = "logrank")
  direct <- analysis_fast(df[df$group %in% c(1, 3), ], control = 1,
                          event.looks = 150)
  p3 <- pw[pw$arm == 3, ]
  p3 <- p3[order(p3$sim), ]
  ok <- direct$reached
  expect_true(all(ok))
  expect_equal(p3$n.event[ok], direct$n.event[ok])
  expect_equal(p3$n.dropout[ok], direct$n.dropout[ok])
  expect_equal(p3$n.pipeline[ok], direct$n.pipeline[ok])
  expect_true(any(p3$n.pipeline > 0))
})

test_that("pairwise_fast: event-driven mode keeps strata columns", {
  df <- make_karm(20)
  df$subgroup <- as.integer(df$accrual_time > 6) + 1L
  expect_error(
    pairwise_fast(df, control = 1, event.looks = 150, primary = 3,
                  stat = "logrank", strata = "subgroup"),
    NA)
})

test_that("pairwise_fast: the adjusted p-value column must be unambiguous", {
  df <- make_karm(20)
  expect_error(
    pairwise_fast(df, control = 1, time.looks = 30, stat = c("logrank", "rmst"),
                  tau = 12, adjust = "bonferroni"),
    "p.col")
  pw <- pairwise_fast(df, control = 1, time.looks = 30,
                      stat = c("logrank", "rmst"), tau = 12,
                      adjust = "bonferroni", p.col = "rmst.p")
  expect_equal(pw$p.adj, pmin(1, 2 * pw$rmst.p))
})

test_that("pairwise_fast: event-driven mode equals re-cutting the data at the primary cutoffs", {
  nsim <- 30
  df <- simdata_fast(nsim = nsim, n = c(120, 120, 120), a.time = c(0, 12),
                     a.rate = 360 / 12, e.median = list(12, 16, 20),
                     d.hazard = 0.02, seed = 919)
  looks <- c(80, 120)
  pw <- pairwise_fast(df, control = 1, event.looks = looks, primary = 3,
                      stat = c("logrank", "rmst"), tau = 12)
  prim <- analysis_fast(df[df$group %in% c(1, 3), ], control = 1,
                        event.looks = looks)
  expect_true(all(prim$reached))
  A <- matrix(prim$cutoff, nrow = nsim, byrow = TRUE)
  big <- max(A) + 1
  for (j in c(2, 3)) {
    for (l in seq_along(looks)) {
      # Reference: administrative censoring at the primary cutoff in R, then
      # a single uncut calendar look.
      sub      <- df[df$group %in% c(1, j), ]
      A_row    <- A[sub$sim, l]
      enrolled <- sub$accrual_time <= A_row
      ended    <- sub$accrual_time + sub$tte <= A_row
      cut_dat       <- sub
      cut_dat$tte   <- pmin(sub$tte, A_row - sub$accrual_time)
      cut_dat$event <- sub$event * as.integer(ended)
      cut_dat       <- cut_dat[enrolled, ]
      ref <- analysis_fast(cut_dat, control = 1, time.looks = big,
                           stat = c("logrank", "rmst"), tau = 12)
      got <- pw[pw$arm == j & pw$look == l, ]
      got <- got[order(got$sim), ]
      expect_equal(got$cutoff, A[, l])
      expect_equal(got$look.value, rep(looks[l], nsim))
      expect_equal(got$logrank.z, ref$logrank.z, tolerance = 1e-10)
      expect_equal(got$rmst.diff, ref$rmst.diff, tolerance = 1e-10)
      expect_equal(got$n.event, ref$n.event)
      expect_equal(got$n.enrolled, ref$n.enrolled)
    }
  }
})

test_that("pairwise_fast: cutoff.looks analyzes every contrast at the supplied cutoffs", {
  df  <- make_karm(40)
  cut <- cutoff_fast(df, event.looks = 150, time.looks = 18)
  pw  <- pairwise_fast(df, control = 1, cutoff.looks = cut, stat = "logrank",
                       side = 1)
  expect_equal(nrow(pw), 2 * 40)
  for (j in 2:3) {
    ref <- analysis_fast(df[df$group %in% c(1, j), ], control = 1,
                         cutoff.looks = cut, stat = "logrank", side = 1)
    got <- pw[pw$arm == j, ]
    got <- got[order(got$sim), ]
    expect_equal(got$cutoff, ref$cutoff)
    expect_equal(got$logrank.z, ref$logrank.z, tolerance = 1e-12)
  }
})

test_that("pairwise_fast: unreached shared cutoffs blank the statistics", {
  df  <- make_karm(10)
  cut <- cutoff_fast(df, event.looks = 100)
  cut[2, 1] <- NA
  pw <- pairwise_fast(df, control = 1, cutoff.looks = cut, stat = "logrank")
  miss <- pw[pw$sim == 2, ]
  expect_false(any(miss$reached))
  expect_true(all(is.na(miss$logrank.z)))
  expect_true(all(is.na(miss$n.event)))
  expect_true(all(!is.na(pw$logrank.z[pw$sim != 2])))
})
