test_that("simdata_fast: structure and schema, no subgroups, two groups", {
  set.seed(1)
  dat <- simdata_fast(
    nsim = 50, n = c(100, 120), a.time = c(0, 12), a.prop = 1,
    e.hazard = list(log(2) / 12, log(2) / 18),
    d.median = list(36, 36), seed = 1
  )

  expect_s3_class(dat, "data.frame")
  expect_equal(nrow(dat), 50 * (100 + 120))
  expect_identical(names(dat),
                   c("sim", "group", "accrual_time", "surv_time",
                     "dropout_time", "tte", "event", "calendar_time"))
  expect_setequal(unique(dat$group), c(1, 2))
  expect_true(all(dat$tte == pmin(dat$surv_time, dat$dropout_time)))
  expect_true(all(dat$event == as.integer(dat$surv_time <= dat$dropout_time)))
  expect_equal(dat$calendar_time, dat$accrual_time + dat$tte, tolerance = 1e-12)
})

test_that("simdata_fast: rows are interleaved in (sim, group) order", {
  set.seed(2)
  nc <- 30; nt <- 40
  dat <- simdata_fast(
    nsim = 10, n = c(nc, nt), a.time = c(0, 10), a.prop = 1,
    e.hazard = list(log(2) / 10, log(2) / 14), seed = 2
  )
  # Within each simulation the first nc rows are control, the next nt treatment.
  for (s in 1:10) {
    block <- dat[dat$sim == s, ]
    expect_equal(block$group, c(rep(1, nc), rep(2, nt)))
  }
})

test_that("simdata_fast: reproducible from seed, sensitive to seed", {
  args <- list(nsim = 20, n = c(50, 50), a.time = c(0, 8), a.prop = 1,
               e.hazard = list(log(2) / 10, log(2) / 12), d.median = list(30, 30))
  a <- do.call(simdata_fast, c(args, list(seed = 123)))
  b <- do.call(simdata_fast, c(args, list(seed = 123)))
  c2 <- do.call(simdata_fast, c(args, list(seed = 124)))

  expect_equal(a$surv_time, b$surv_time, tolerance = 0)
  expect_equal(a$accrual_time, b$accrual_time, tolerance = 0)
  expect_false(isTRUE(all.equal(a$surv_time, c2$surv_time)))
})

test_that("simdata_fast: single-group mode", {
  set.seed(3)
  dat <- simdata_fast(
    nsim = 30, n = 80, a.time = c(0, 12), a.prop = 1,
    e.hazard = log(2) / 12, seed = 3
  )
  expect_equal(nrow(dat), 30 * 80)
  expect_setequal(unique(dat$group), 1)
  expect_true(all(is.infinite(dat$dropout_time)))  # no dropout supplied
  expect_true(all(dat$event == 1L))
})

test_that("simdata_fast: exponential survival recovers the hazard", {
  set.seed(4)
  haz <- log(2) / 15
  dat <- simdata_fast(
    nsim = 1, n = 50000, a.time = c(0, 1), a.prop = 1,
    e.hazard = haz, seed = 4
  )
  # MLE of an exponential rate with full follow-up is 1 / mean(time).
  expect_equal(1 / mean(dat$surv_time), haz, tolerance = 0.03)
})

test_that("simdata_fast: piecewise survival changes hazard at breakpoints", {
  set.seed(5)
  dat <- simdata_fast(
    nsim = 1, n = 80000, a.time = c(0, 1), a.prop = 1,
    e.hazard = c(0.05, 0.30), e.time = c(0, 6, Inf), seed = 5
  )
  # Early hazard low, late hazard high: empirical hazard over [0,6) should be
  # well below that over [6, Inf).
  early <- mean(dat$surv_time < 6)
  expect_gt(mean(dat$surv_time >= 6), 0.2)
  # Crude check: fraction failing before 6 is approximately 1 - exp(-0.05 * 6).
  expect_equal(early, 1 - exp(-0.05 * 6), tolerance = 0.03)
})

test_that("simdata_fast: prevalence is recovered (single factor)", {
  set.seed(6)
  prev <- c(0.6, 0.4)
  dat <- simdata_fast(
    nsim = 1, n = 60000, a.time = c(0, 1), a.prop = 1,
    e.hazard = log(2) / 12, prevalence = prev, seed = 6
  )
  expect_true("subgroup" %in% names(dat))
  tab <- prop.table(table(dat$subgroup))
  expect_equal(as.numeric(tab), prev, tolerance = 0.02)
})

test_that("simdata_fast: independent multi-factor prevalence", {
  set.seed(7)
  prev <- list(c(0.5, 0.5), c(0.7, 0.3))
  dat <- simdata_fast(
    nsim = 1, n = 60000, a.time = c(0, 1), a.prop = 1,
    e.hazard = log(2) / 12, prevalence = prev, seed = 7
  )
  expect_true(all(c("subgroup1", "subgroup2") %in% names(dat)))
  expect_equal(as.numeric(prop.table(table(dat$subgroup1))),
               c(0.5, 0.5), tolerance = 0.02)
  expect_equal(as.numeric(prop.table(table(dat$subgroup2))),
               c(0.7, 0.3), tolerance = 0.02)
})

test_that("simdata_fast: per-cell subgroup hazards are recovered", {
  set.seed(8)
  # One factor, two cells, with cell-specific survival hazards.
  dat <- simdata_fast(
    nsim = 1, n = 80000, a.time = c(0, 1), a.prop = 1,
    e.hazard = list(0.05, 0.20),       # per-cell hazards
    prevalence = c(0.5, 0.5), seed = 8
  )
  h1 <- 1 / mean(dat$surv_time[dat$subgroup == 1])
  h2 <- 1 / mean(dat$surv_time[dat$subgroup == 2])
  expect_equal(h1, 0.05, tolerance = 0.01)
  expect_equal(h2, 0.20, tolerance = 0.02)
})

test_that("simdata_fast: two-group subgroup simultaneous recovery", {
  set.seed(9)
  dat <- simdata_fast(
    nsim = 1, n = c(40000, 40000), a.time = c(0, 1), a.prop = 1,
    e.hazard = list(log(2) / 12, log(2) / 20),
    prevalence = c(0.6, 0.4), seed = 9
  )
  expect_setequal(unique(dat$group), c(1, 2))
  expect_true("subgroup" %in% names(dat))
  # Group hazards recovered marginally over subgroups.
  hc <- 1 / mean(dat$surv_time[dat$group == 1])
  ht <- 1 / mean(dat$surv_time[dat$group == 2])
  expect_equal(hc, log(2) / 12, tolerance = 0.02)
  expect_equal(ht, log(2) / 20, tolerance = 0.02)
  # Prevalence recovered within each group.
  expect_equal(as.numeric(prop.table(table(dat$subgroup[dat$group == 1]))),
               c(0.6, 0.4), tolerance = 0.02)
})

test_that("simdata_fast: dropout produces censoring", {
  set.seed(10)
  no_drop <- simdata_fast(nsim = 1, n = 20000, a.time = c(0, 1), a.prop = 1,
                          e.hazard = log(2) / 12, seed = 10)
  with_drop <- simdata_fast(nsim = 1, n = 20000, a.time = c(0, 1), a.prop = 1,
                            e.hazard = log(2) / 12, d.median = 12, seed = 10)
  expect_true(all(no_drop$event == 1L))
  expect_lt(mean(with_drop$event), 1)
  expect_gt(mean(with_drop$event), 0.3)
})

test_that("simdata_fast: fixed.alloc gives deterministic subgroup counts", {
  set.seed(11)
  n_one <- 100
  dat <- simdata_fast(
    nsim = 5, n = n_one, a.time = c(0, 1), a.prop = 1,
    e.hazard = log(2) / 12, prevalence = c(0.6, 0.4),
    fixed.alloc = TRUE, seed = 11
  )
  # Each simulation has identical, deterministic cell counts.
  for (s in 1:5) {
    cnt <- as.numeric(table(dat$subgroup[dat$sim == s]))
    expect_equal(cnt, c(60, 40))
  }
})

test_that("simdata_fast: degenerate prevalence equals no-subgroup output", {
  a <- simdata_fast(nsim = 10, n = c(50, 50), a.time = c(0, 8), a.prop = 1,
                    e.hazard = list(log(2) / 10, log(2) / 12), seed = 77)
  b <- simdata_fast(nsim = 10, n = c(50, 50), a.time = c(0, 8), a.prop = 1,
                    e.hazard = list(log(2) / 10, log(2) / 12),
                    prevalence = c(1), seed = 77)
  expect_equal(a$surv_time, b$surv_time, tolerance = 0)
  expect_equal(a$tte, b$tte, tolerance = 0)
})

test_that("simdata_fast: two-group with per-cell control hazards (nested list)", {
  set.seed(12)
  # Control has cell-specific hazards (three cells), treatment shares one.
  dat <- simdata_fast(
    nsim = 1, n = c(60000, 60000), a.time = c(0, 1), a.prop = 1,
    e.hazard = list(list(0.10, 0.08, 0.06), 0.05),
    prevalence = c(0.4, 0.3, 0.3), seed = 12
  )
  expect_setequal(unique(dat$group), c(1, 2))
  # Control cell hazards recovered.
  for (cl in 1:3) {
    h <- 1 / mean(dat$surv_time[dat$group == 1 & dat$subgroup == cl])
    expect_equal(h, c(0.10, 0.08, 0.06)[cl], tolerance = 0.02,
                 info = paste("control cell", cl))
  }
  # Treatment hazard is common across cells.
  ht <- 1 / mean(dat$surv_time[dat$group == 2])
  expect_equal(ht, 0.05, tolerance = 0.01)
})

test_that("simdata_fast: input validation", {
  expect_error(
    simdata_fast(nsim = 5, n = c(1, 2, 3), a.time = c(0, 1), a.prop = 1,
                 e.hazard = log(2) / 12),
    "scalar .total N., a vector of length 2")
  expect_error(
    simdata_fast(nsim = 5, n = 50, a.time = c(0, 1), a.rate = 1),
    "One of 'e.hazard' or 'e.median'")
  expect_error(
    simdata_fast(nsim = 5, n = 50, a.time = c(0, 1), a.prop = 1,
                 e.hazard = log(2) / 12, e.median = 12),
    "exactly one of")
  expect_error(
    simdata_fast(nsim = 5, n = 50, a.time = c(0, 1), a.rate = 5,
                 e.hazard = log(2) / 12),
    "a.rate")
  expect_error(
    simdata_fast(nsim = 5, n = 50, a.time = c(0, 1), a.prop = 1,
                 e.hazard = c(0.05, 0.2)),
    "e.time")
})

# ------------------------------------------------------------------ #
#  Accrual specification: absolute rate, final-interval completion,
#  proportions, and deterministic per-interval counts
# ------------------------------------------------------------------ #

test_that("simdata_fast: a.rate places deterministic per-interval counts", {
  dat <- simdata_fast(
    nsim = 7, n = 360, a.time = c(0, 12, 24), a.rate = c(10, 20),
    e.median = 18, seed = 1)
  per_sim <- split(dat$accrual_time, dat$sim)
  cnt1 <- vapply(per_sim, function(a) sum(a < 12), integer(1))
  expect_true(all(cnt1 == 120))
  expect_true(all(vapply(per_sim, length, integer(1)) == 360))
  expect_true(all(dat$accrual_time >= 0 & dat$accrual_time < 24))
})

test_that("simdata_fast: a.rate completes the final interval from the total", {
  dat <- simdata_fast(
    nsim = 5, n = 500, a.time = c(0, 12), a.rate = c(20, 30),
    e.median = 18, seed = 1)
  per_sim <- split(dat$accrual_time, dat$sim)
  cnt1 <- vapply(per_sim, function(a) sum(a < 12), integer(1))
  expect_true(all(cnt1 == 240))
  end_time <- 12 + 260 / 30
  expect_lt(max(dat$accrual_time), end_time + 1e-8)
  expect_gt(max(dat$accrual_time), 12)
})

test_that("simdata_fast: a.prop distributes the total by proportion", {
  dat <- simdata_fast(
    nsim = 4, n = 100, a.time = c(0, 6, 12), a.prop = c(0.3, 0.7),
    e.median = 18, seed = 1)
  per_sim <- split(dat$accrual_time, dat$sim)
  cnt1 <- vapply(per_sim, function(a) sum(a < 6), integer(1))
  expect_true(all(cnt1 == 30))
  expect_true(all(dat$accrual_time <= 12))
})

test_that("simdata_fast: accrual counts split across two groups", {
  dat <- simdata_fast(
    nsim = 3, n = c(180, 180), a.time = c(0, 12, 24), a.rate = c(10, 20),
    e.median = list(18, 18), seed = 2)
  per_sim <- split(dat$accrual_time, dat$sim)
  cnt1 <- vapply(per_sim, function(a) sum(a < 12), integer(1))
  expect_true(all(cnt1 == 120))
})

test_that("simdata_fast: accrual specification validation", {
  # absolute rate inconsistent with the total
  expect_error(
    simdata_fast(nsim = 2, n = 500, a.time = c(0, 12, 24), a.rate = c(10, 20),
                 e.median = 18),
    "implies")
  # a.rate with an invalid length
  expect_error(
    simdata_fast(nsim = 2, n = 50, a.time = c(0, 1), a.rate = c(1, 2, 3),
                 e.median = 18),
    "a.rate")
  # exactly one of a.rate / a.prop must be supplied
  expect_error(
    simdata_fast(nsim = 2, n = 50, a.time = c(0, 12), a.rate = 50 / 12,
                 a.prop = 1, e.median = 18),
    "exactly one of 'a.rate'")
  expect_error(
    simdata_fast(nsim = 2, n = 50, a.time = c(0, 12), e.median = 18),
    "exactly one of 'a.rate'")
  # a.prop length must match the number of intervals
  expect_error(
    simdata_fast(nsim = 2, n = 50, a.time = c(0, 6, 12), a.prop = 1,
                 e.median = 18),
    "a.prop")
})

test_that("simdata_fast: a scalar n is split without losing subjects", {
  df <- simdata_fast(nsim = 2, n = 7, a.time = c(0, 7), a.rate = 1,
                     e.hazard = list(0.1, 0.08), seed = 1)
  expect_equal(nrow(df), 14)
  expect_equal(as.vector(table(df$group[df$sim == 1])), c(4, 3))
  df2 <- simdata_fast(nsim = 2, n = 301, a.time = c(0, 12), a.prop = 1,
                      e.hazard = list(0.1, 0.08), seed = 1)
  expect_equal(nrow(df2), 602)
  expect_error(simdata_fast(nsim = 2, n = 100, alloc = c(1, 1, 1),
                            a.time = c(0, 12), a.prop = 1,
                            e.hazard = list(0.1, 0.08)), "alloc")
  expect_error(simdata_fast(nsim = 2, n = 100, alloc = c(1, -1),
                            a.time = c(0, 12), a.prop = 1,
                            e.hazard = list(0.1, 0.08)), "alloc")
})

test_that("simdata_fast: hazards and piecewise breakpoints are validated", {
  expect_error(simdata_fast(nsim = 2, n = 50, a.time = c(0, 1), a.prop = 1,
                            e.hazard = c(0.05, 0.2), e.time = c(3, 6, Inf)),
               "start at 0")
  expect_error(simdata_fast(nsim = 2, n = 50, a.time = c(0, 1), a.prop = 1,
                            e.hazard = c(0.05, 0.2, 0.1),
                            e.time = c(0, 12, 6, Inf)),
               "strictly increasing")
  expect_error(simdata_fast(nsim = 2, n = 50, a.time = c(0, 1), a.prop = 1,
                            e.hazard = -0.1), "non-negative")
  expect_error(simdata_fast(nsim = 2, n = 50, a.time = c(0, 1), a.prop = 1,
                            e.median = -5), "non-negative")
})

test_that("simdata_fast: fixed.alloc subgroups are not tied to the accrual order", {
  df <- simdata_fast(nsim = 200, n = 150, a.time = c(0, 6, 12),
                     a.rate = c(5, 20), e.hazard = 0.05,
                     prevalence = c(0.2, 0.8), fixed.alloc = TRUE, seed = 11)
  # Fixed counts: 30 and 120 subjects in every simulation.
  cnt <- tapply(df$subgroup == 1L, df$sim, sum)
  expect_true(all(cnt == 30))
  # The 30 subjects enrolled in [0, 6] are a random draw from the 150, so on
  # average 20% of them are in subgroup 1 (it was 100% before the fix).
  early <- df$accrual_time < 6
  expect_lt(abs(mean(df$subgroup[early] == 1L) - 0.2), 0.03)
})

test_that("simdata_fast: fixed.alloc counts are robust to floating-point shares", {
  # Hand calculation: 100 * 0.29 is 28.999... in floating point but the
  # intended count is 29.
  df <- simdata_fast(nsim = 1, n = 100, a.time = c(0, 12), a.prop = 1,
                     e.hazard = 0.05, prevalence = c(0.71, 0.29),
                     fixed.alloc = TRUE, seed = 12)
  expect_equal(as.vector(table(df$subgroup)), c(71, 29))
})

test_that("simdata_fast: scalar n with alloc equals per-group n (issue #2)", {
  # A scalar total split by 'alloc' must reproduce the per-group
  # specification exactly, including subgroups with per-cell hazards.
  args <- list(nsim = 20, a.time = c(0, 12), a.rate = 300 / 12,
               e.hazard = list(list(0.10, 0.08, 0.06), 0.05),
               prevalence = c(0.5, 0.3, 0.2), seed = 3)
  ref <- do.call(simdata_fast, c(list(n = c(150, 150)), args))
  expect_identical(do.call(simdata_fast,
                           c(list(n = 300, alloc = c(1, 1)), args)), ref)
  # A nested (per-group, per-cell) list also signals two groups without alloc.
  expect_identical(do.call(simdata_fast, c(list(n = 300), args)), ref)
  # Unequal allocation.
  ref2 <- do.call(simdata_fast, c(list(n = c(100, 200)), args))
  expect_identical(do.call(simdata_fast,
                           c(list(n = 300, alloc = c(1, 2)), args)), ref2)
})

test_that("simdata_fast: with subgroups, alloc decides per-group versus per-cell lists", {
  base <- list(nsim = 5, a.time = c(0, 12), a.prop = 1,
               e.hazard = list(0.10, 0.05), prevalence = c(0.5, 0.5),
               seed = 21)
  # Scalar n without alloc: one group, per-cell hazards (unchanged behavior).
  one <- do.call(simdata_fast, c(list(n = 200), base))
  expect_equal(unique(one$group), 1L)
  # Scalar n with alloc: two groups, per-group hazards, as with n = c(100, 100).
  two <- do.call(simdata_fast, c(list(n = 200, alloc = c(1, 1)), base))
  ref <- do.call(simdata_fast, c(list(n = c(100, 100)), base))
  expect_identical(two, ref)
})

test_that("simdata_fast: a shared (non-list) dropout works with two groups and subgroups", {
  base <- list(nsim = 5, n = c(100, 100), a.time = c(0, 12), a.prop = 1,
               e.hazard = list(0.10, 0.05), prevalence = c(0.5, 0.5),
               seed = 22)
  shared <- do.call(simdata_fast, c(base, list(d.hazard = 0.01)))
  per_group <- do.call(simdata_fast, c(base, list(d.hazard = list(0.01, 0.01))))
  expect_identical(shared, per_group)
})

test_that("simdata_fast: list lengths are checked instead of silently truncated", {
  expect_error(
    simdata_fast(nsim = 2, n = c(50, 50), a.time = c(0, 12), a.prop = 1,
                 e.hazard = list(0.10, 0.08, 0.05)),
    "list of length 2")
  expect_error(
    simdata_fast(nsim = 2, n = 100, a.time = c(0, 12), a.prop = 1,
                 e.hazard = list(0.10, 0.08, 0.05), prevalence = c(0.5, 0.5)),
    "one element per subgroup cell")
  expect_error(
    simdata_fast(nsim = 2, n = c(50, 50), a.time = c(0, 12), a.prop = 1,
                 e.hazard = list(list(0.10, 0.08), 0.05),
                 prevalence = c(0.3, 0.3, 0.4)),
    "one element per subgroup cell")
})

test_that("simdata_fast: per-cell breakpoints are used in a one-group simulation", {
  # Cell 2 has hazard 0.2 on [0, 3) and 0.1 afterwards, so
  # P(T > 6) = exp(-0.2 * 3 - 0.1 * 3) = exp(-0.9), about 0.407 (hand
  # calculation). Using cell 1's breakpoint 6 instead would give exp(-1.2),
  # about 0.301.
  df <- simdata_fast(nsim = 1, n = 40000, a.time = c(0, 1), a.prop = 1,
                     e.hazard = list(c(0.1, 0.05), c(0.2, 0.1)),
                     e.time = list(c(0, 6, Inf), c(0, 3, Inf)),
                     prevalence = c(0.5, 0.5), seed = 23)
  expect_equal(unique(df$group), 1L)
  p2 <- mean(df$surv_time[df$subgroup == 2L] > 6)
  expect_lt(abs(p2 - exp(-0.9)), 0.02)
})

test_that("simdata_fast: group-specific prevalence must have matching factor levels", {
  expect_error(
    simdata_fast(nsim = 2, n = c(50, 50), a.time = c(0, 12), a.prop = 1,
                 e.hazard = 0.1,
                 prevalence = list(control   = list(c(0.5, 0.5), c(0.3, 0.3, 0.4)),
                                   treatment = list(c(0.3, 0.3, 0.4), c(0.5, 0.5)))),
    "levels per factor")
})

test_that("simdata_fast: stream gives reproducible and distinct batches", {
  args <- list(nsim = 5, n = c(50, 50), a.time = c(0, 6), a.rate = 100 / 6,
               e.median = list(10, 12), seed = 77)
  b1  <- do.call(simdata_fast, c(args, list(stream = 1)))
  b2  <- do.call(simdata_fast, c(args, list(stream = 2)))
  b1b <- do.call(simdata_fast, c(args, list(stream = 1)))
  # Same stream, same data, whatever was generated in between.
  expect_identical(b1b, b1)
  expect_false(isTRUE(all.equal(b1$surv_time, b2$surv_time)))
  # Without 'stream' the seeding of the versions before 'stream' existed is
  # unchanged: the first rows below were printed by the mrct-regional-
  # consistency vignette of FastSurvival 1.0.0 on CRAN (seed = 2; the accrual
  # times are drawn for all simulated trials before the survival times, so
  # nsim must be 10000).
  old <- simdata_fast(nsim = 10000, n = c(25, 25), a.time = c(3, 12.5),
                      a.prop = 1, e.median = list(4.3, 5.811), seed = 2)
  expect_equal(head(old$accrual_time, 6),
               c(8.314255, 10.349985, 5.726526, 7.179403, 5.562820, 9.503599),
               tolerance = 1e-6)
  expect_equal(head(old$surv_time, 6),
               c(1.589398, 1.604216, 3.582047, 3.952740, 1.644677, 15.446380),
               tolerance = 1e-6)
  # The illness-death and multi-arm modes use the same stream mechanism.
  id_args <- list(nsim = 3, n = c(40, 40), a.time = c(0, 6), a.rate = 80 / 6,
                  h01.hazard = list(0.10, 0.07), h02.hazard = list(0.05, 0.04),
                  seed = 5)
  expect_identical(do.call(simdata_fast, c(id_args, list(stream = 3))),
                   do.call(simdata_fast, c(id_args, list(stream = 3))))
  k_args <- list(nsim = 3, n = c(40, 40, 40), a.time = c(0, 6),
                 a.rate = 120 / 6, e.median = list(10, 12, 14), seed = 6)
  k1 <- do.call(simdata_fast, c(k_args, list(stream = 1)))
  k2 <- do.call(simdata_fast, c(k_args, list(stream = 2)))
  expect_false(isTRUE(all.equal(k1$surv_time, k2$surv_time)))

  expect_error(simdata_fast(nsim = 2, n = 20, a.time = c(0, 6),
                            a.rate = 20 / 6, e.median = 10, stream = 1),
               "requires 'seed'")
  expect_error(do.call(simdata_fast, c(args, list(stream = -1))),
               "non-negative whole number")
  expect_error(do.call(simdata_fast, c(args, list(stream = 1.5))),
               "non-negative whole number")
})

test_that("simdata_fast: a subject without a finite event or dropout time is censored", {
  # Cure model: zero hazard after month 24 and no dropout.
  df <- simdata_fast(nsim = 20, n = c(100, 100), a.time = c(0, 12),
                     a.rate = 200 / 12, e.hazard = list(c(0.1, 0), c(0.1, 0)),
                     e.time = c(0, 24, Inf), seed = 5)
  never <- !is.finite(df$surv_time)
  expect_true(any(never))
  expect_true(all(df$surv_time[!never] <= 24))
  expect_lt(abs(mean(never) - exp(-0.1 * 24)), 0.02)
  expect_equal(df$event, as.integer(!never))
  expect_true(all(df$tte[never] == Inf))

  # Event-driven looks count only finite events, and an unreached look keeps
  # the never-event subjects in follow-up (not counted as dropouts).
  n_ev <- tapply(df$event, df$sim, sum)
  looks <- c(min(n_ev), max(n_ev) + 1)
  res <- analysis_fast(df, control = 1, event.looks = looks)
  l1 <- res$look == 1
  expect_true(all(res$reached[l1]) && all(is.finite(res$cutoff[l1])))
  expect_true(all(!res$reached[!l1]) && all(is.na(res$cutoff[!l1])))
  expect_true(all(res$n.dropout[!l1] == 0))
  cut <- cutoff_fast(df, event.looks = looks)
  expect_equal(unname(cut[, 1]), res$cutoff[l1])
  expect_true(all(is.na(cut[, 2])))
})

test_that("simdata_fast: a single-level prevalence gives its subgroup column", {
  args <- list(nsim = 3, n = c(20, 20), a.time = c(0, 6), a.rate = 40 / 6,
               e.median = list(12, 18), seed = 3)
  a <- do.call(simdata_fast, c(args, list(prevalence = 1)))
  b <- do.call(simdata_fast, args)
  expect_true("subgroup" %in% names(a))
  expect_true(all(a$subgroup == 1L))
  # A single cell consumes no random numbers.
  expect_identical(a[names(b)], b)
  res <- analysis_fast(a, control = 1, time.looks = 12, by.subgroup = TRUE)
  expect_setequal(unique(res$population), c("overall", "subgroup_1"))
})

test_that("simdata_fast: 'n' must contain positive whole numbers", {
  base <- list(nsim = 2, a.time = c(0, 6), a.rate = 40 / 6,
               e.median = list(12, 18))
  for (bad in list(c(20, 20.5), c(20, 0), c(-20, 60), c(20, NA), 0, 10.5,
                   "20")) {
    expect_error(do.call(simdata_fast, c(base, list(n = bad))),
                 "positive whole numbers")
  }
})

test_that("simdata_fast: stream = 0 gives the same data as no stream", {
  args <- list(nsim = 3, n = c(30, 30), a.time = c(0, 6), a.rate = 10,
               e.median = list(12, 18), seed = 11)
  expect_identical(do.call(simdata_fast, c(args, list(stream = 0))),
                   do.call(simdata_fast, args))
  expect_false(identical(do.call(simdata_fast, c(args, list(stream = 1))),
                         do.call(simdata_fast, args)))
})

test_that("simdata_fast: a whole number up to rounding error is not truncated", {
  # 90 * 0.7 is 62.99999999999999 in floating point.
  n <- 90 * c(0.3, 0.7)
  expect_false(n[2] == 63)
  d <- simdata_fast(nsim = 2, n = n, a.time = c(0, 12), a.prop = 1,
                    e.median = list(12, 18), seed = 1)
  tab <- table(d$sim, d$group)
  expect_equal(as.vector(tab[, "1"]), c(27L, 27L))
  expect_equal(as.vector(tab[, "2"]), c(63L, 63L))
  # With fixed subgroup sizes, every simulated trial has the fixed cell counts
  # of the whole number (9 and 18 of 27, 21 and 42 of 63).
  df <- simdata_fast(nsim = 3, n = n, a.time = c(0, 12), a.prop = 1,
                     e.median = list(12, 18), prevalence = c(1, 2) / 3,
                     fixed.alloc = TRUE, seed = 1)
  for (s in 1:3) {
    ds <- df[df$sim == s, ]
    expect_equal(as.vector(table(ds$subgroup[ds$group == 1L])), c(9L, 18L))
    expect_equal(as.vector(table(ds$subgroup[ds$group == 2L])), c(21L, 42L))
  }
})

test_that("simdata_fast: error messages name the piecewise argument", {
  expect_error(
    simdata_fast(nsim = 1, n = c(10, 10), a.time = c(0, 1), a.prop = 1,
                 e.median = list(10, 10), d.hazard = c(0.01, 0.02)),
    "'d.time' must be supplied")
  expect_error(
    simdata_fast(nsim = 1, n = c(10, 10), a.time = c(0, 1), a.prop = 1,
                 e.median = list(10, 10), d.hazard = c(0.01, 0.02),
                 d.time = c(0, 6)),
    "length\\(d.time\\)")
  expect_error(
    simdata_fast(nsim = 1, n = c(10, 10), a.time = c(0, 1), a.prop = 1,
                 h01.hazard = c(0.1, 0.2), h02.hazard = 0.05),
    "'h01.time' must be supplied")
})

test_that("simdata_fast: nsim must be a whole number", {
  expect_error(simdata_fast(nsim = 2.5, n = c(10, 10), a.time = c(0, 1),
                            a.prop = 1, e.median = list(10, 10)),
               "'nsim' must be a positive whole number")
})

test_that("simdata_fast rounds nsim and stream within rounding error of a whole number", {
  ns <- 100 * 0.29  # 28.999999999999996, truncated to 28 before 1.2.0
  d <- simdata_fast(nsim = ns, n = c(5, 5), a.time = c(0, 1), a.prop = 1,
                    e.median = list(10, 10), seed = 1)
  expect_equal(length(unique(d$sim)), 29L)
  args <- list(nsim = 2, n = c(5, 5), a.time = c(0, 1), a.prop = 1,
               e.median = list(10, 10), seed = 1)
  s_f  <- do.call(simdata_fast, c(args, list(stream = ns)))
  s_29 <- do.call(simdata_fast, c(args, list(stream = 29)))
  s_28 <- do.call(simdata_fast, c(args, list(stream = 28)))
  expect_identical(s_f, s_29)
  expect_false(identical(s_f, s_28))
})
