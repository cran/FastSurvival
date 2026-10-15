make_print_data <- function() {
  simdata_fast(nsim = 1, n = c(60, 60), a.time = c(0, 6), a.prop = 1,
               e.median = list(8, 14), d.hazard = 0.02, seed = 501)
}

test_that("print.milestone_fast labels the p-value by the side and method", {
  d <- make_print_data()
  x2 <- milestone_fast(d$tte, d$event, d$group, control = 1, tau = 6)
  expect_output(v <- withVisible(print(x2)), "Pr(>|z|)", fixed = TRUE)
  expect_false(v$visible)
  expect_identical(v$value, x2)
  # One-sided: the upper tail for "wald" and "mover", the lower tail for
  # "loglog", whose statistic is negative when treatment is better.
  for (m in c("wald", "mover", "loglog")) {
    x1 <- milestone_fast(d$tte, d$event, d$group, control = 1, tau = 6,
                         side = 1, method = m)
    out1 <- utils::capture.output(print(x1))
    lab <- if (m == "loglog") "Pr(<z)" else "Pr(>z)"
    expect_true(any(grepl(lab, out1, fixed = TRUE)), info = m)
    expect_false(any(grepl("Pr(>|z|)", out1, fixed = TRUE)), info = m)
  }
})
