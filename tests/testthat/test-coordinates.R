context("Manipulating coordinates")

test_that("makeFriendlyCoord works on strings and numerics", {
    expect_equal(makeFriendlyCoord("chr1", 1000, 2000), "chr1:1000-2000")
    expect_equal(makeFriendlyCoord("chr1", "1000", "2000"), "chr1:1000-2000")
})
