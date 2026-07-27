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
               "exactly one of 'event.looks' or 'time.looks'")
  expect_error(pairwise_fast(df, control = 1, event.looks = 100, time.looks = 30),
               "exactly one of 'event.looks' or 'time.looks'")
  expect_error(pairwise_fast(df, control = 1, event.looks = 100),
               "'primary' must name")
  expect_error(pairwise_fast(df, control = 9, time.looks = 30),
               "not among the group labels")
  expect_error(pairwise_fast(df, control = 1, arms = c(1, 2), time.looks = 30),
               "must not include the control")
  expect_error(pairwise_fast(df, control = 1, time.looks = 30, by.subgroup = TRUE),
               "does not support 'by.subgroup'")
})
