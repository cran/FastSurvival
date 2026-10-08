make_scn <- function() {
  gen_scenario_fast(
    scenarios = list(
      "Proportional"   = list(e.median = list(12, 18)),
      "Delayed effect" = list(
        e.hazard = list(log(2) / 12, c(log(2) / 12, log(2) / 22)),
        e.time   = c(0, 6, Inf)
      ),
      "Crossing"       = list(
        e.hazard = list(log(2) / 12, c(log(2) / 7, log(2) / 24)),
        e.time   = c(0, 5, Inf),
        null     = TRUE
      )
    ),
    shared = list(n = c(150, 150), a.time = c(0, 12), a.rate = 300 / 12)
  )
}

test_that("gen_scenario_fast merges shared and scenario arguments", {
  scn <- make_scn()
  expect_s3_class(scn, "scenario_fast")
  expect_equal(names(scn$scenarios),
               c("Proportional", "Delayed effect", "Crossing"))
  a1 <- scn$scenarios[[1]]$args
  expect_equal(a1$n, c(150, 150))
  expect_equal(a1$e.median, list(12, 18))
  # Scenario arguments take precedence over shared ones.
  scn2 <- gen_scenario_fast(list(list(e.median = 10, n = 80)),
                            shared = list(n = 200, a.time = c(0, 6),
                                          a.rate = 80 / 6))
  expect_equal(scn2$scenarios[[1]]$args$n, 80)
  # 'null' and 'label' are interpreted and removed from the arguments.
  expect_true(scn$scenarios[[3]]$null)
  expect_false(scn$scenarios[[1]]$null)
  expect_false("null" %in% names(scn$scenarios[[3]]$args))
  expect_equal(scn$shared$n, c(150, 150))
})

test_that("gen_scenario_fast resolves labels in the documented order", {
  sc <- list(a = list(e.median = 10, label = "from field"),
             b = list(e.median = 12),
             list(e.median = 14))
  scn <- gen_scenario_fast(sc)
  expect_equal(names(scn$scenarios), c("from field", "b", "Scenario 3"))
  expect_false("label" %in% names(scn$scenarios[[1]]$args))
  scn_l <- gen_scenario_fast(sc, labels = c("L1", "L2", "L3"))
  # A 'label' field still overrides 'labels'.
  expect_equal(names(scn_l$scenarios), c("from field", "L2", "L3"))
})

test_that("gen_scenario_fast arguments reproduce a direct simdata_fast call", {
  scn <- make_scn()
  s2  <- scn$scenarios[["Delayed effect"]]$args
  d1  <- do.call(simdata_fast, c(s2, list(nsim = 3, seed = 9)))
  d2  <- simdata_fast(nsim = 3, n = c(150, 150), a.time = c(0, 12),
                      a.rate = 300 / 12,
                      e.hazard = list(log(2) / 12,
                                      c(log(2) / 12, log(2) / 22)),
                      e.time = c(0, 6, Inf), seed = 9)
  expect_identical(d1, d2)
})

test_that("gen_scenario_fast analytic helpers match hand calculations", {
  # Piecewise hazard 0.1 on [0, 6) and 0.05 afterward:
  # H(3) = 0.3, H(10) = 0.6 + 0.2 = 0.8.
  expect_equal(pwe_cumhaz(c(3, 10), c(0.1, 0.05), c(0, 6, Inf)), c(0.3, 0.8))
  expect_equal(pwe_hazard(c(3, 10), c(0.1, 0.05), c(0, 6, Inf)), c(0.1, 0.05))
  expect_equal(pwe_cumhaz(4, 0.2, NULL), 0.8)
  # Median of an exponential curve with median 12, read by interpolation.
  tt <- seq(0, 48, length.out = 4801)
  expect_equal(curve_median(tt, exp(-log(2) / 12 * tt)), 12, tolerance = 1e-4)
  expect_true(is.na(curve_median(tt, rep(0.9, length(tt)))))
  # Hazard ratio of the delayed-effect scenario: 1 before month 6 and
  # 12 / 22 afterward.
  ev <- scenario_eval(make_scn()$scenarios[["Delayed effect"]]$args, c(3, 10))
  expect_equal(ev$hr, c(1, 12 / 22))
  expect_true(ev$two)
})

test_that("gen_scenario_fast validates its input", {
  expect_error(gen_scenario_fast(list()), "non-empty list")
  expect_error(gen_scenario_fast("a"), "non-empty list")
  expect_error(gen_scenario_fast(list(1)), "must be a named list")
  expect_error(gen_scenario_fast(list(list(n = 10))), "no 'e.hazard'")
  expect_error(gen_scenario_fast(list(list(e.median = 10)), shared = 1),
               "'shared'")
  expect_error(gen_scenario_fast(list(list(e.median = 10)),
                                 labels = c("a", "b")),
               "length\\(labels\\)")
})
