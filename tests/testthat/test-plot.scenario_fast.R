test_that("plot.scenario_fast draws the scenarios and restores the layout", {
  scn <- gen_scenario_fast(
    scenarios = list(
      "Proportional" = list(e.median = list(12, 18)),
      "Crossing"     = list(
        e.hazard = list(log(2) / 12, c(log(2) / 7, log(2) / 24)),
        e.time   = c(0, 5, Inf)
      ),
      "Single"       = list(e.median = 10)
    ),
    shared = list(n = c(150, 150), a.time = c(0, 12), a.rate = 300 / 12)
  )
  grDevices::pdf(tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  mf0 <- graphics::par("mfrow")
  v <- withVisible(plot(scn))
  expect_false(v$visible)
  expect_identical(v$value, scn)
  expect_equal(graphics::par("mfrow"), mf0)
  expect_error(plot(scn, which = "Crossing", tmax = 30, hr_max = 2,
                    legend_pos = "bottomleft"), NA)
  expect_error(plot(scn, which = c(1, 3), mfrow = c(1, 2)), NA)
  expect_error(plot(scn, which = integer(0)), "no scenarios to plot")
  bad <- gen_scenario_fast(list(list(e.hazard = list(0.1, c(0.1, 0.05)))))
  expect_error(plot(bad), "'e.time' is required")
})
