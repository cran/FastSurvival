make_print_data <- function() {
  simdata_fast(nsim = 1, n = c(60, 60), a.time = c(0, 6), a.prop = 1,
               e.median = list(8, 14), d.hazard = 0.02, seed = 501)
}

test_that("print.survfit_fast prints the estimate and returns x invisibly", {
  d <- make_print_data()
  o <- order(d$tte)
  x <- survfit_fast(d$tte[o], d$event[o], t_eval = 6, conf.type = "log-log")
  expect_output(v <- withVisible(print(x)),
                "Kaplan-Meier survival estimate \\(single time point\\)")
  expect_false(v$visible)
  expect_identical(v$value, x)
  expect_output(print(x), "t = 6")
  expect_output(print(x), "Confidence interval type: log-log")
  xn <- x
  xn[] <- NA_real_
  expect_output(print(xn), "Estimate not available")
})
