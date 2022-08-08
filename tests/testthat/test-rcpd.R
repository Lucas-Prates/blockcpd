test_that("rcpd sampling", {
  error_msg = "Input error! Change point vector entries must vary from 1 to ncol-1"
  expect_error(rcpd(nrow = 1, ncol = 5, changepoints = c(1, 10)))

  changepoints = c(2, 5, 8)
  td = rcpd(nrow = 10, ncol = 10, changepoints = changepoints, family = "normal",
            parameter = list(mean = c(1, 2, 3, 4), var = c(0.1, 0.2 , 0.3, 0.4)))
  expect_equal(td$changepoints, changepoints)

  warning_msg = "NA"
  expect_warning(rcpd(parameters = list(p = c(0.1, 2))), warning_msg)
})
