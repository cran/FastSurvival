make_print_data <- function() {
  simdata_fast(nsim = 1, n = c(60, 60), a.time = c(0, 6), a.prop = 1,
               e.median = list(8, 14), d.hazard = 0.02, seed = 501)
}

test_that("print.maxcombo_fast prints the components and the combined test", {
  d <- make_print_data()
  x <- maxcombo_fast(d$tte, d$event, d$group, control = 1, side = 1,
                     rho = c(0, 0), gamma = c(0, 1))
  expect_output(v <- withVisible(print(x)),
                "Max-combo weighted log-rank test \\(two-group\\)")
  expect_false(v$visible)
  expect_identical(v$value, x)
  expect_output(print(x), "N = 120,  control = 1", fixed = TRUE)
  expect_output(print(x), "FH(0,1)", fixed = TRUE)
  expect_output(print(x), "(one-sided)", fixed = TRUE)
  set.seed(1)
  x2 <- maxcombo_fast(d$tte, d$event, d$group, control = 1, side = 2,
                      rho = c(0, 0), gamma = c(0, 1))
  expect_output(print(x2), "(two-sided)", fixed = TRUE)
  xn <- x
  xn["statistic"] <- NA_real_
  expect_output(print(xn), "Test statistic not available")
})
