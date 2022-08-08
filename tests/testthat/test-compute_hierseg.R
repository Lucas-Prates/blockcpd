test_that("hierseg c(1, 1, 1, 1)", {
  model = fit_blockcpd(t(as.matrix(c(1, 1, 1, 1))),
                       lambda = 0, method = "hierseg")
  expect_equal(model$ncp, 0)
  expect_equal(model$neg_loglike, 0)
})

test_that("hierseg c(0, 1, 1, 0)", {
  model = fit_blockcpd(c(0, 1, 1, 0),
                       lambda = 0, method = "hierseg")
  expect_equal(model$changepoints, c(1, 3))
  expect_equal(model$ncp, 2)
  expect_equal(model$neg_loglike, 0)
})

test_that("hierseg c(1, 2, 3, 4) normal", {
  model = fit_blockcpd(c(1, 2, 3, 4), family = "normal",
                       lambda = 0, method = "hierseg")
  expect_equal(model$changepoints, c(2))
  expect_equal(model$ncp, 1)
})
