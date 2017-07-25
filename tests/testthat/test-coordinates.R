context("Manipulating coordinates")

test_that("Conversion between strings and GRanges objects works", {

    string <- "chr1:1000-2000"
    gr <- GRanges(seqnames = "chr1",
                  ranges = IRanges(start = 1000, end = 2000))

    expect_equal(makeBrowserCoord("chr1", 1000, 2000), string)
    expect_equal(makeBrowserCoord("chr1", "1000", "2000"), string)

    expect_equal(coordToGRanges(string), gr)
    expect_equal(GRangesToCoord(gr), string)

})
