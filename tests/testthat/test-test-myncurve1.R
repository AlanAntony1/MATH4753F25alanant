test_that("mu is correct", {
  result <- myncurve(0,1,1)
  expect_equal(result$mu, 0)
})

test_that("sigma is correct", {
  result <- myncurve(0,1,1)
  expect_equal(result$sigma, 1)
})

test_that("probability is correct", {
  result <- myncurve(0,1,1)
  expected_prob <- pnorm(1, mean = 0, sd = 1)
  expect_equal(result$probability, expected_prob)
})
