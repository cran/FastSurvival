test_that("print.scenario_fast reports the analytic medians and hazard ratios", {
  scn <- gen_scenario_fast(
    scenarios = list(
      "Proportional"   = list(e.median = list(12, 18)),
      "Delayed effect" = list(
        e.hazard = list(log(2) / 12, c(log(2) / 12, log(2) / 22)),
        e.time   = c(0, 6, Inf)
      ),
      "Crossing"       = list(
        e.hazard = list(log(2) / 12, c(log(2) / 7, log(2) / 24)),
        e.time   = c(0, 5, Inf)
      )
    ),
    shared = list(n = c(150, 150), a.time = c(0, 12), a.rate = 300 / 12)
  )
  expect_output(v <- withVisible(print(scn)),
                "A scenario_fast object with 3 scenarios")
  expect_false(v$visible)
  expect_identical(v$value, scn)

  tab <- do.call(rbind, lapply(scn$scenarios, scenario_summary_row,
                               tmax = 48))
  expect_equal(tab$N, rep(300L, 3))
  # Hand calculations: medians 12 and 18 (proportional); the delayed-effect
  # treatment median solves 6 / 12 + (t - 6) / 22 = 1, so t = 17; the crossing
  # treatment median solves 5 / 7 + (t - 5) / 24 = 1, so t = 5 + 48 / 7.
  expect_equal(tab$Median_C, rep(12, 3), tolerance = 2e-3)
  expect_equal(tab$Median_T, c(18, 17, 5 + 48 / 7), tolerance = 2e-3)
  expect_equal(tab$HR_start, round(c(12 / 18, 1, 12 / 7), 3))
  expect_equal(tab$HR_end, round(c(12 / 18, 12 / 22, 12 / 24), 3))
  expect_equal(tab$Crossing, c(FALSE, FALSE, TRUE))

  # A single-group scenario has no treatment median or hazard ratio.
  one <- gen_scenario_fast(list(One = list(e.median = 10, n = 100)))
  expect_output(print(one), "1 scenario ")
  row <- scenario_summary_row(one$scenarios[[1]], tmax = 40)
  expect_true(is.na(row$Median_T))
  expect_true(is.na(row$HR_start))
  expect_true(is.na(row$Crossing))
})
