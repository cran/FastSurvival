make_print_data <- function() {
  simdata_fast(nsim = 1, n = c(60, 60), a.time = c(0, 6), a.prop = 1,
               e.median = list(8, 14), d.hazard = 0.02, seed = 501)
}

test_that("print.rmw_fast prints the components and the combined test", {
  d <- make_print_data()
  x <- rmw_fast(d$tte, d$event, d$group, control = 1, side = 2)
  expect_output(v <- withVisible(print(x)),
                "Robust modestly-weighted log-rank test \\(two-group\\)")
  expect_false(v$visible)
  expect_identical(v$value, x)
  expect_output(print(x), "N = 120,  control = 1,  s_star = 0.5", fixed = TRUE)
  expect_output(print(x), "Null correlation", fixed = TRUE)
  expect_output(print(x), "two-sided p-value", fixed = TRUE)
  x1 <- rmw_fast(d$tte, d$event, d$group, control = 1, side = 1)
  expect_output(print(x1), "one-sided p-value", fixed = TRUE)
  xn <- x
  xn[1] <- NA_real_
  expect_output(print(xn), "Test statistic not available")
})
