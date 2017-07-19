context("Local processing functions")

test_that("retrievePeaks finds peaks in the query region", {

    metadata <- data.frame(Sample = c("A", "B"),
                           TestMark = c("../test_data/basic.bed",
                                        "../test_data/basic.bed"),
                           stringsAsFactors = FALSE)

    pks <- list(A = GRanges(seqnames = rep("chr1", 3),
                            ranges = IRanges(start = c(100, 150, 500),
                                             end = c(200, 250, 600))),
                B = GRanges(seqnames = rep("chr1", 3),
                            ranges = IRanges(start = c(100, 150, 500),
                                             end = c(200, 250, 600))))

    pks_in_region <- list(A = GRanges(seqnames = rep("chr1", 2),
                                      ranges = IRanges(start = c(100, 150),
                                                       end = c(200, 250))),
                          B = GRanges(seqnames = rep("chr1", 2),
                                      ranges = IRanges(start = c(100, 150),
                                                       end = c(200, 250))))

    pks_in_region$A <- `elementMetadata<-`(x = pks_in_region$A,
                        value = data.frame(length = as.integer(c(101, 101))))
    pks_in_region$B <- `elementMetadata<-`(x = pks_in_region$B,
                        value = data.frame(length = as.integer(c(101, 101))))

    region <- GRanges(seqnames = "chr1",
                      ranges = IRanges(start = 100, end = 300))

    expect_equal(retrievePeaks(pks, metadata, region),
                 localPeaks(region, pks_in_region, c("A", "B")))

    region2 <- GRanges(seqnames = "chr2",
                      ranges = IRanges(start = 100, end = 300))

    # No peaks
    expect_warning(retrievePeaks(pks, metadata, region2))

})
