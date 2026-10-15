test_that("plot.kmcurve_fast draws the number-at-risk table and a band", {
  set.seed(98)
  n  <- 60
  t0 <- rexp(n, log(2) / 12)
  t1 <- rexp(n, log(2) / 18)
  cc <- runif(2 * n, 0, 30)
  time  <- pmin(c(t0, t1), cc)
  event <- as.integer(c(t0, t1) <= cc)
  group <- rep(c(1, 2), each = n)
  fit2 <- kmcurve_fast(time, event, group, control = 1)
  fit1 <- kmcurve_fast(time[group == 1], event[group == 1])

  grDevices::pdf(tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  mar0 <- graphics::par("mar")
  mf0  <- graphics::par("mfrow")

  v <- withVisible(plot(fit2, risk_table = TRUE))
  expect_false(v$visible)
  expect_identical(v$value, fit2)
  # The margins and the layout are restored
  expect_equal(graphics::par("mar"), mar0)
  expect_equal(graphics::par("mfrow"), mf0)

  expect_error(plot(fit2, hr = TRUE, risk_table = TRUE,
                    risk_times = c(0, 10, 20, 100), conf.style = "band",
                    group_labels = c("Placebo", "Active drug")), NA)
  expect_error(plot(fit2, risk_table = TRUE,
                    group_labels = factor(c("Placebo", "Drug"))), NA)
  expect_error(plot(fit1, risk_table = TRUE, conf.style = "band",
                    conf.type = "log-log"), NA)
  expect_error(plot(fit2, conf.int = FALSE, conf.style = "band"), NA)
  expect_error(plot(fit2, risk_table = TRUE, base_size = 2, hr = TRUE), NA)
  expect_equal(graphics::par("mar"), mar0)

  expect_error(plot(fit2, conf.style = "shade"), "should be one of")
  expect_error(plot(fit2, risk_table = TRUE, risk_times = -1),
               "non-negative")
  expect_error(plot(fit2, risk_table = TRUE, risk_times = c(1, NA)),
               "non-negative")
})
