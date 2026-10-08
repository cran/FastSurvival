make_print_data <- function() {
  simdata_fast(nsim = 1, n = c(60, 60), a.time = c(0, 6), a.prop = 1,
               e.median = list(8, 14), d.hazard = 0.02, seed = 501)
}

test_that("print.ahsw_fast prints the average hazards and contrasts", {
  d <- make_print_data()
  x <- ahsw_fast(d$tte, d$event, d$group, control = 1, tau = 10)
  expect_output(v <- withVisible(print(x)),
                "Average hazard with survival weight \\(two-group\\)")
  expect_false(v$visible)
  expect_identical(v$value, x)
  expect_output(print(x), "tau = 10,  control = 1", fixed = TRUE)
  expect_output(print(x), "ratio (treatment / control)", fixed = TRUE)
  expect_output(print(x), "difference (treatment - control)", fixed = TRUE)
  x1 <- ahsw_fast(d$tte, d$event, d$group, control = 1, tau = 10, side = 1)
  expect_output(print(x1), "alternative = one.sided")
  xn <- x
  xn[["ah.ctrl"]] <- NA_real_
  expect_output(print(xn), "Estimate not available")
})
