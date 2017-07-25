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

    uniq_pks <- GRanges(seqnames = rep("chr1", 2),
                        ranges = IRanges(start = c(100, 130),
                                         end = c(150, 180)))

    expect_equal(getUniquePeaks(
        GRanges(seqnames = rep("chr1", 4),
                ranges = IRanges(start = c(100, 110, 120, 130),
                end = c(150, 160, 170, 180))), 0.5), uniq_pks)
})


test_that("getSamplePeakProfile correctly models peak presence", {

    pks <- GRanges(seqnames = rep("chr1", 4),
                    ranges = IRanges(start = c(100, 110, 120, 130),
                                     end = c(150, 160, 170, 180)))

    pks2 <- GRanges(seqnames = rep("chr4", 4),
                   ranges = IRanges(start = c(100, 110, 120, 130),
                                    end = c(150, 160, 170, 180)))

    uniq_pks <- GRanges(seqnames = rep("chr1", 2),
                        ranges = IRanges(start = c(100, 130),
                                         end = c(150, 180)))

    expect_equal(getSamplePeakProfile(pks, uniq_pks, 0.9),
                 data.frame(X1 = TRUE, X2 = TRUE))

    expect_equal(suppressWarnings(getSamplePeakProfile(pks2, uniq_pks, 0.9)),
                 data.frame(X1 = FALSE, X2 = FALSE))

})


test_that("position-aware feature matrix construction works", {

    a <- GRanges(seqnames = rep("chr1", 2),
                ranges = IRanges(start = c(10, 2050), end = c(20, 2500)))
    b <- GRanges(seqnames = rep("chr1", 2),
                ranges = IRanges(start = c(10, 2000), end = c(22, 2600)))
    c <- GRanges(seqnames = rep("chr1", 2),
                ranges = IRanges(start = c(10, 2020), end = c(1000, 2450)))
    d <- GRanges(seqnames = rep("chr1", 2),
                ranges = IRanges(start = c(10, 2030), end = c(999, 2700)))

    pks <- list(a = a, b = b, c = c, d = d)

    lp <- localPeaks(region = GRanges(seqnames = "chr1",
                                    ranges = IRanges(start = 1, 3000)),
                    peaks = pks,
                    samples = c("a", "b", "c", "d"))

    position_matrix <- data.frame(ft1 = c(TRUE, TRUE, FALSE, FALSE),
                                ft2 = c(TRUE, TRUE, TRUE, TRUE),
                                ft3 = c(FALSE, FALSE, TRUE, TRUE))

    colnames(position_matrix) = c("chr1:10-20",
                                    "chr1:2050-2500",
                                    "chr1:10-1000")
    rownames(position_matrix) = lpkSamples(lp)
    attr(position_matrix, "features") <- getUniquePeaks(Reduce("c", pks), 0.5)

    expect_equal(positionAware(lp, 0.5), position_matrix)

})
