test_that("step_vertices traces a right-continuous step function", {
  # Hand calculation: 1 on [0, 1), 0.8 on [1, 3), and 0.5 at 3
  v <- step_vertices(c(0, 1, 3), c(1, 0.8, 0.5))
  expect_equal(v$x, c(0, 1, 1, 3, 3))
  expect_equal(v$y, c(1, 1, 0.8, 0.8, 0.5))
  # A single point
  v1 <- step_vertices(2, 0.4)
  expect_equal(v1$x, 2)
  expect_equal(v1$y, 0.4)
})
