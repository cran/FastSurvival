# Analysis output of six simulated trials at two looks
pss_data <- function() {
  df <- data.frame(
    sim       = rep(1:6, each = 2),
    look      = rep(1:2, 6),
    logrank.z = c(-2.9, -3.1, -1.0, -2.2, 0.5, -1.0,
                  -2.0, -2.5, -0.2, 0.3, -1.5, -2.1),
    cox.coef  = c(-0.4, -0.5, -0.1, -0.3, 0.3, 0.1,
                  -0.2, -0.3, 0.25, 0.2, -0.2, -0.3),
    n.event   = rep(c(50, 100), 6),
    cutoff    = rep(c(12, 24), 6)
  )
  df$logrank.p <- pnorm(df$logrank.z)
  df
}

test_that("plot.simsummary_fast draws the stopping probabilities", {
  df  <- pss_data()
  res <- simsummary_fast(df, eff.col = "logrank.z", efficacy = c(-2.8, -1.96),
                         fut.col = "cox.coef", futility = c(0.2, NA))
  dp   <- rbind(cbind(df, population = "A"), cbind(df, population = "B"))
  resp <- simsummary_fast(dp, eff.col = "logrank.z", efficacy = c(-2.8, -1.96),
                          fut.col = "cox.coef", futility = c(0.2, NA))
  resP <- simsummary_fast(df, p.col = "logrank.p", alpha = c(0.0026, 0.025))

  grDevices::pdf(tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  mf0 <- graphics::par("mfrow")

  v <- withVisible(plot(res))
  expect_false(v$visible)
  expect_identical(v$value, res)
  expect_error(plot(res, xaxis = "events"), NA)
  expect_error(plot(res, what = "stage", xaxis = "time"), NA)
  expect_error(plot(resp), NA)
  expect_error(plot(resp, what = "stage", col = c("red", "blue")), NA)
  expect_equal(graphics::par("mfrow"), mf0)
  expect_error(plot(resP, what = "stage"), NA)
  expect_error(plot(resP, ylim = c(0, 1), main = "OC"), NA)
})

test_that("plot.simsummary_fast rejects incomplete results", {
  df  <- pss_data()
  res <- simsummary_fast(df, eff.col = "logrank.z", efficacy = c(-2.8, -1.96))
  expect_error(plot(res[res$look != "1", ]), "keep all their looks")
  expect_error(plot(res, what = "bars"), "should be one of")
  res_nc <- simsummary_fast(df[, setdiff(names(df), "cutoff")],
                            eff.col = "logrank.z", efficacy = c(-2.8, -1.96))
  expect_error(plot(res_nc, xaxis = "time"), "cutoff.mean")
  df_na <- df
  df_na$cutoff <- NA_real_
  res_na <- simsummary_fast(df_na, eff.col = "logrank.z",
                            efficacy = c(-2.8, -1.96))
  expect_error(plot(res_na, xaxis = "time"), "finite values")
  # A result of FastSurvival 1.1.0 stores no looks or rows
  old <- res
  bd <- attr(old, "boundary")
  bd$looks <- NULL
  bd$rows <- NULL
  attr(old, "boundary") <- bd
  expect_error(plot(old), "earlier version of FastSurvival")
})
