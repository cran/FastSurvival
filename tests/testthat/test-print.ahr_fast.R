make_print_data <- function() {
  simdata_fast(nsim = 1, n = c(60, 60), a.time = c(0, 6), a.prop = 1,
               e.median = list(8, 14), d.hazard = 0.02, seed = 501)
}

test_that("print.ahr_fast prints the hazard shares and the average hazard ratio", {
  d <- make_print_data()
  x <- ahr_fast(d$tte, d$event, d$group, control = 1, tau = 10)
  expect_output(v <- withVisible(print(x)),
                "Kalbfleisch-Prentice average hazard ratio \\(two-group\\)")
  expect_false(v$visible)
  expect_identical(v$value, x)
  expect_output(print(x), "average hazard ratio (treatment / control)",
                fixed = TRUE)
  expect_output(print(x), "(log scale: z =", fixed = TRUE)
  x1 <- ahr_fast(d$tte, d$event, d$group, control = 1, tau = 10, side = 1)
  expect_output(print(x1), "alternative = one.sided")
  xn <- x
  xn$ahr <- NA_real_
  expect_output(print(xn), "Estimate not available")
})

test_that("print.ahr_fast labels the p-value by the side of the test", {
  d <- make_print_data()
  x2 <- ahr_fast(d$tte, d$event, d$group, control = 1, tau = 10)
  expect_output(print(x2), "Pr(>|z|)", fixed = TRUE)
  x1 <- ahr_fast(d$tte, d$event, d$group, control = 1, tau = 10, side = 1)
  out1 <- utils::capture.output(print(x1))
  expect_true(any(grepl("Pr(<z)", out1, fixed = TRUE)))
  expect_false(any(grepl("Pr(>|z|)", out1, fixed = TRUE)))
})
