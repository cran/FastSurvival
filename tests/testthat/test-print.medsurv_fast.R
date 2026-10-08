make_print_data <- function() {
  simdata_fast(nsim = 1, n = c(60, 60), a.time = c(0, 6), a.prop = 1,
               e.median = list(8, 14), d.hazard = 0.02, seed = 501)
}

test_that("print.medsurv_fast prints the two-group and single-group summaries", {
  d <- make_print_data()
  x <- medsurv_fast(d$tte, d$event, d$group, control = 1)
  expect_output(v <- withVisible(print(x)),
                "Median survival time \\(two-group\\)")
  expect_false(v$visible)
  expect_identical(v$value, x)
  expect_output(print(x), "method = km,  alternative = two.sided", fixed = TRUE)
  expect_output(print(x), "difference (treatment - control)", fixed = TRUE)
  xn <- x
  xn["median.control"] <- NA_real_
  expect_output(print(xn), "Estimate not available")

  s <- medsurv_fast(d$tte, d$event, method = "nph")
  expect_output(print(s), "Median survival time \\(single-group\\)")
  expect_output(print(s), "method = nph", fixed = TRUE)
  sn <- s
  sn["median"] <- NA_real_
  expect_output(print(sn), "Estimate not available")
})
