test_that("Compare identical change-point sets",{
    metrics = compare_model(c(1, 2, 3), c(1, 2, 3), 10)
    expect_equal(metrics$haus, 0)
    expect_equal(metrics$rand, 1)
    expect_equal(metrics$symdiff, 0)
    expect_equal(metrics$jaccard, 1)
})

test_that("Raised error on missing ncol argument",{
  error_msg =  "Error! No blockcpd models were passed and ncol was not provided!"
  expect_error(compare_model(c(1, 2, 3), c(1, 2, 3)), error_msg)
})

test_that("Compare different change-points sets",{
  metrics = compare_model(c(1, 5), c(1, 2, 3), 10)
  expect_equal(metrics$haus, 2)
  expect_equal(round(metrics$rand, 3), round(30/45, 3))
  expect_equal(metrics$symdiff, 3)
  expect_equal(metrics$jaccard, 0.4)
})
