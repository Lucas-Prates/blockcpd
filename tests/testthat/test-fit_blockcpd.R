test_that("fit_blockcpd bernoulli", {
  model = fit_blockcpd(c(0, 1, 0, 1), lambda = Inf)
  expect_equal(model$ncp, 0)

  model = fit_blockcpd(c(0, 1, 0, 1), lambda = 0)
  expect_equal(model$ncp, 3)
})

test_that("fit_blockcpd normal", {
  model = fit_blockcpd(c(1, 2, 3, 4, 5, 6), family = "normal", method = "dynseg",
                       lambda = 0)
  expect_equal(model$changepoints, c(2, 4))

  model = fit_blockcpd(c(1, 2, 3, 4, 5, 6), family = "normal", method = "dynseg",
                       min_block_size = 3, lambda = 0)
  expect_equal(model$changepoints, c(3))
})


test_that("fit_blockcpd poisson", {
  model = fit_blockcpd(c(1, 2, 3, 4),
                       lambda = 0, family = "poisson")
  expect_equal(model$ncp, 3)

  model = fit_blockcpd(c(1, 2, 3, 4),
                       lambda = 0, family = "poisson",
                       min_block_size = 2)
  expect_equal(model$ncp, 1)
  expect_equal(model$changepoints, c(2))
})

test_that("fit_blockcpd exponential", {
  model = fit_blockcpd(c(1, 2, 3, 4),
                       lambda = 0, family = "exponential")
  expect_equal(model$ncp, 3)

  model = fit_blockcpd(c(1, 2, 3, 4),
                       lambda = 0, family = "exponential",
                       min_block_size = 2)
  expect_equal(model$ncp, 1)
  expect_equal(model$changepoints, c(2))
})

test_that("fit_blockcpd input check", {
  error_msg = "Input error! The 'lambda' argument must be non-negative"
  expect_error(fit_blockcpd(c(1, 1), lambda = -1), error_msg)

  error_msg = "Input error! The 'min_block_size' argument ranges from 1 to ncol!"
  expect_error(fit_blockcpd(c(1, 1), min_block_size = 0), error_msg)

  error_msg = "Input error! The data passed has only one column!"
  expect_error(fit_blockcpd(1), error_msg)

  error_msg = "Input error! The 'family' argument provided is not implemented!"
  expect_error(fit_blockcpd(c(1, 1), family = "banana"), error_msg)
})
