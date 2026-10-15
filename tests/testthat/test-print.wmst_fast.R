make_print_data <- function() {
  simdata_fast(nsim = 1, n = c(60, 60), a.time = c(0, 6), a.prop = 1,
               e.median = list(8, 14), d.hazard = 0.02, seed = 501)
}

test_that("print.wmst_fast prints the two-group and single-group summaries", {
  d <- make_print_data()
  x <- wmst_fast(d$tte, d$event, d$group, control = 1, tau1 = 2, tau2 = 10)
  expect_output(v <- withVisible(print(x)),
                "Window mean survival time \\(two-group\\)")
  expect_false(v$visible)
  expect_identical(v$value, x)
  expect_output(print(x), "window = [2, 10],  control = 1", fixed = TRUE)
  expect_output(print(x), "difference (treatment - control)", fixed = TRUE)
  s <- wmst_fast(d$tte, d$event, tau1 = 2, tau2 = 10, side = 1)
  expect_output(print(s), "Window mean survival time \\(single-group\\)")
  expect_output(print(s), "window = [2, 10]", fixed = TRUE)
})

test_that("print.wmst_fast labels the one-sided p-value Pr(>z)", {
  d <- make_print_data()
  x2 <- wmst_fast(d$tte, d$event, d$group, control = 1, tau1 = 2, tau2 = 10)
  expect_output(print(x2), "Pr(>|z|)", fixed = TRUE)
  x1 <- wmst_fast(d$tte, d$event, d$group, control = 1, tau1 = 2, tau2 = 10,
                  side = 1)
  out1 <- utils::capture.output(print(x1))
  expect_true(any(grepl("Pr(>z)", out1, fixed = TRUE)))
  expect_false(any(grepl("Pr(>|z|)", out1, fixed = TRUE)))
})
