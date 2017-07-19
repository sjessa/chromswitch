context("Manipulating coordinates")

test_that("makeBrowserCoord works on strings and numerics", {
    expect_equal(makeBrowserCoord("chr1", 1000, 2000), "chr1:1000-2000")
    expect_equal(makeBrowserCoord("chr1", "1000", "2000"), "chr1:1000-2000")
})
