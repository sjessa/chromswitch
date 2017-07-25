context("Position-aware strategy for feature matrix construction")

test_that("pReciprocalOverlap properly distinguishes regions", {

    a <- GRanges(seqnames = "chr11",
                 ranges = IRanges(start = 112829468, end = 112834468))
    b <- GRanges(seqnames = "chr11",
                 ranges = IRanges(start = 112829468, end = 113834468))
    c <- GRanges(seqnames = "chr11",
                 ranges = IRanges(start = 112829968, end = 112834968))
    d <- GRanges(seqnames = "chr11",
                 ranges = IRanges(start = 112830468, end = 112833468))
    e <- GRanges(seqnames = "chr12",
                 ranges = IRanges(start = 112830468, end = 112833468))

    expect_equal(pReciprocalOverlap(a, b, 0.9), FALSE)
    expect_equal(pReciprocalOverlap(b, a, 0.9), FALSE)
    expect_equal(pReciprocalOverlap(a, c, 0.9), TRUE)
    expect_equal(pReciprocalOverlap(a, c, 0.95), FALSE)
    expect_equal(pReciprocalOverlap(a, c, 0.5), TRUE)
    expect_equal(pReciprocalOverlap(a, d, 0.9), FALSE)
    expect_warning(pReciprocalOverlap(a, e, 0.9),
                   "no sequence levels in common")
    expect_equal(suppressWarnings(pReciprocalOverlap(a, e, 0.9)), FALSE)

})


test_that("getUniquePeaks collapses a list of peaks by p-reciprocal overlap", {

    a <- GRanges(seqnames = "chr11",
                 ranges = IRanges(start = 112829468, end = 112834468))
    b <- GRanges(seqnames = "chr11",
                 ranges = IRanges(start = 112829468, end = 113834468))
    c <- GRanges(seqnames = "chr11",
                 ranges = IRanges(start = 112829968, end = 112834968))
    d <- GRanges(seqnames = "chr11",
                 ranges = IRanges(start = 112830468, end = 112833468))
    e <- GRanges(seqnames = "chr12",
                 ranges = IRanges(start = 112830468, end = 112833468))

    pks <- suppressWarnings(c(a, b, c, d, e))

    expect_equal(getUniquePeaks(pks, 0.9), suppressWarnings(c(a, b, d, e)))

    pks2 <- GRanges(seqnames = rep("chr1", 4),
                    ranges = IRanges(start = c(100, 110, 120, 130),
                                     end = c(150, 160, 170, 180)))

    expect_equal(getUniquePeaks(
        GRanges(seqnames = rep("chr1", 4),
                ranges = IRanges(start = c(100, 110, 120, 130),
                end = c(150, 160, 170, 180))), 0.5),
        GRanges(seqnames = rep("chr1", 2),
                ranges = IRanges(start = c(100, 130),
                                end = c(150, 180))))
})
