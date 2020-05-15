test_that("computeS works", {
    expect_error(computeS(1:5), "matrix")
    expect_equal(computeS(matrix(rep(1, 5), 5, 1)),
                 matrix(1/5, 5, 5))
})
