make_print_data <- function() {
  simdata_fast(nsim = 1, n = c(60, 60), a.time = c(0, 6), a.prop = 1,
               e.median = list(8, 14), d.hazard = 0.02, seed = 501)
}

test_that("print.coxph_fast labels the p-value by the side of the test", {
  d <- make_print_data()
  x2 <- coxph_fast(d$tte, d$event, d$group, control = 1)
  expect_output(v <- withVisible(print(x2)), "Pr(>|z|)", fixed = TRUE)
  expect_false(v$visible)
  expect_identical(v$value, x2)
  x1 <- coxph_fast(d$tte, d$event, d$group, control = 1, side = 1)
  out1 <- utils::capture.output(print(x1))
  expect_true(any(grepl("Pr(<z)", out1, fixed = TRUE)))
  expect_false(any(grepl("Pr(>|z|)", out1, fixed = TRUE)))
  expect_true(any(grepl("one.sided", out1, fixed = TRUE)))
})
