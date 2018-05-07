context("test cal.eigen.fit")

test_that("test cal.eigen.fit",{

  X <- c(2.7668157, 1.8931580, 1.8859049, 1.0575971, 1.0329016, 0.9095669,
         0.6395215, 0.5073142, 0.3205904, 0.3125665)

  res <- cal.eigen.fit(X)

  expected_vec <- c(0.379450739, 0.003838584, 0.578408288, 0.023627539,
                    0.127158683, 0.352248264, 0.231589695, 0.458966294,
                    0.025346914)

  expect_length(res,9)
  expect_equal(round(res,5), round(expected_vec,5))
  expect_type(res, "double")

})
