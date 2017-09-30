context("localPeaks class")

test_that("localPeaks access works", {

    pks_in_region <- list(A = GRanges(seqnames = rep("chr1", 2),
                                      ranges = IRanges(start = c(100, 150),
                                                       end = c(200, 250))),
                          B = GRanges(seqnames = rep("chr1", 2),
                                      ranges = IRanges(start = c(100, 150),
                                                       end = c(200, 250))))

    region <- GRanges(seqnames = "chr1",
                      ranges = IRanges(start = 100, end = 300))

    lpk <- localPeaks(region, pks_in_region, c("A", "B"))

    expect_equal(region(lpk), region)
    expect_equal(samples(lpk), c("A", "B"))
    expect_equal(peaks(lpk), pks_in_region)

})


test_that("isEmpty finds empty localPeaks", {

    region <- GRanges(seqnames = "chr1",
                      ranges = IRanges(start = 100, end = 300))
    pks    <- GRanges(seqnames = c(), ranges = IRanges(start = c(), end = c()))
    pks_in_region <- list(A = GRanges(seqnames = rep("chr1", 2),
                                      ranges = IRanges(start = c(100, 150),
                                                       end = c(200, 250))),
                          B = GRanges(seqnames = rep("chr1", 2),
                                      ranges = IRanges(start = c(100, 150),
                                                       end = c(200, 250))))
    expect_equal(isEmpty(localPeaks(region,
                                     list(A = pks, B = pks),
                                     c("A", "B"))),
                 TRUE)

    expect_equal(isEmpty(localPeaks(region,
                                     list(A = pks_in_region, B = pks_in_region),
                                     c("A", "B"))),
                 FALSE)

})
