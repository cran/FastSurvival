make_print_data <- function() {
  simdata_fast(nsim = 1, n = c(60, 60), a.time = c(0, 6), a.prop = 1,
               e.median = list(8, 14), d.hazard = 0.02, seed = 501)
}

test_that("print.rmst_fast prints the two-group and single-group summaries", {
  d <- make_print_data()
  x <- rmst_fast(d$tte, d$event, d$group, control = 1, tau = 10)
  expect_output(v <- withVisible(print(x)),
                "Restricted mean survival time \\(two-group\\)")
  expect_false(v$visible)
  expect_identical(v$value, x)
  expect_output(print(x), "tau = 10,  control = 1", fixed = TRUE)
  expect_output(print(x), "difference (treatment - control)", fixed = TRUE)
  expect_output(print(x), "ratio (treatment / control)", fixed = TRUE)
  x1 <- rmst_fast(d$tte, d$event, d$group, control = 1, tau = 10, side = 1)
  expect_output(print(x1), "alternative = one.sided")

  s <- rmst_fast(d$tte, d$event, tau = 10)
  expect_output(print(s), "Restricted mean survival time \\(single-group\\)")
  sn <- s
  sn[["rmst"]] <- NA_real_
  expect_output(print(sn), "Estimate not available")
})

test_that("print.rmst_fast labels the p-value by the side of the test", {
  d <- make_print_data()
  x2 <- rmst_fast(d$tte, d$event, d$group, control = 1, tau = 10)
  expect_output(print(x2), "Pr(>|z|)", fixed = TRUE)
  x1 <- rmst_fast(d$tte, d$event, d$group, control = 1, tau = 10, side = 1)
  out1 <- utils::capture.output(print(x1))
  expect_true(any(grepl("Pr(>z)", out1, fixed = TRUE)))
  expect_false(any(grepl("Pr(>|z|)", out1, fixed = TRUE)))
})
