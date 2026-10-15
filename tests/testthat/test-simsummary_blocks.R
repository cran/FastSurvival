# Analysis output of six simulated trials at two looks
ssb_data <- function() {
  data.frame(
    sim       = rep(1:6, each = 2),
    look      = rep(1:2, 6),
    logrank.z = c(-2.9, -3.1, -1.0, -2.2, 0.5, -1.0,
                  -2.0, -2.5, -0.2, 0.3, -1.5, -2.1),
    cox.coef  = c(-0.4, -0.5, -0.1, -0.3, 0.3, 0.1,
                  -0.2, -0.3, 0.25, 0.2, -0.2, -0.3),
    n.event   = rep(c(50, 100), 6),
    cutoff    = rep(c(12, 24), 6)
  )
}

test_that("simsummary_blocks returns the look rows of each block", {
  df  <- ssb_data()
  res <- simsummary_fast(df, eff.col = "logrank.z", efficacy = c(-2.8, -1.96),
                         fut.col = "cox.coef", futility = c(0.2, NA))
  b <- simsummary_blocks(res)
  expect_length(b, 1L)
  expect_identical(attr(b[[1L]], "label"), "")
  expect_identical(b[[1L]]$look, c("1", "2"))
  expect_equal(b[[1L]]$cum.reject, res$cum.reject[1:2])

  # Populations and arms are separate blocks, labeled as in the print method
  dp   <- rbind(cbind(df, population = "A"), cbind(df, population = "B"))
  resp <- simsummary_fast(dp, eff.col = "logrank.z", efficacy = c(-2.8, -1.96))
  bp   <- simsummary_blocks(resp)
  expect_identical(vapply(bp, attr, "", "label"),
                   c("Population: A", "Population: B"))
  da   <- rbind(cbind(df, arm = 2), cbind(df, arm = 3))
  resa <- simsummary_fast(da, eff.col = "logrank.z", efficacy = c(-2.8, -1.96))
  expect_identical(vapply(simsummary_blocks(resa), attr, "", "label"),
                   c("Arm: 2", "Arm: 3"))
  # One population of several is a valid result of its own
  expect_length(simsummary_blocks(resp[resp$population == "B", ]), 1L)
})

test_that("simsummary_blocks rejects selections that the report rejects", {
  df  <- ssb_data()
  res <- simsummary_fast(df, eff.col = "logrank.z", efficacy = c(-2.8, -1.96))
  expect_null(simsummary_blocks(res[res$look != "2", ]))
  expect_null(simsummary_blocks(res[c(2, 1, 3), ]))
  expect_null(simsummary_blocks(res[, c("population", "look", "cum.reject")]))
  # Rows of another result combined by rbind()
  other <- simsummary_fast(df, eff.col = "logrank.z", efficacy = c(-2, -1))
  other$population <- "B"
  expect_null(simsummary_blocks(rbind(res, other)))
  # The overall row is not needed
  expect_length(simsummary_blocks(res[res$look != "overall", ]), 1L)
})
