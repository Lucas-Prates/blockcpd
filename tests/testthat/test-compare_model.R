test_that("Compare identical change-point sets",
          {
            expect_equal(compare_model(c(1, 2, 3), c(1, 2, 3), 10),
                         list(haus = 0,
                              rand = 1,
                              symdiff = 0,
                              jaccard = 1)
            )
          }
)

test_that("Raised error on missing ncol argument",
          {
            expect_error(compare_model(c(1, 2, 3), c(1, 2, 3)),
                         "Error! No blockcpd models were passed and ncol was not provided!",
            )
          }
)

test_that("Compare different change-points sets",
          {
            expect_equal(round(unlist(compare_model(c(1, 5), c(1, 2, 3), 10)), 3),
                         round(unlist(list(haus = 2,
                                           rand = 30/45,
                                           symdiff = 3,
                                           jaccard = 0.4)), 3)
            )
          }
)
