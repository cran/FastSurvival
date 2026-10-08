make_print_data <- function() {
  simdata_fast(nsim = 1, n = c(60, 60), a.time = c(0, 6), a.prop = 1,
               e.median = list(8, 14), d.hazard = 0.02, seed = 501)
}

test_that("print.wkm_fast prints the weighted Kaplan-Meier test", {
  d <- make_print_data()
  x <- wkm_fast(d$tte, d$event, d$group, control = 1, weight = "sqrtPF",
                side = 1)
  expect_output(v <- withVisible(print(x)),
                "Weighted Kaplan-Meier test \\(Pepe-Fleming, two-group\\)")
  expect_false(v$visible)
  expect_identical(v$value, x)
  expect_output(print(x), "weight = sqrtPF,  alternative = one.sided",
                fixed = TRUE)
  expect_output(print(x), "weighted difference (treatment - control)",
                fixed = TRUE)
  xn <- x
  xn["z"] <- NA_real_
  expect_output(print(xn), "Test statistic not available")
})
