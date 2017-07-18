context("Reading BED files")

test_that("loadSampleBed handles BED files properly", {

    expect_equal(loadSampleBed(bed_file = "../test_data/basic.bed"),
                 GRanges(seqnames = rep("chr1", 3),
                         ranges = IRanges(start = c(100, 150, 500),
                                          end = c(200, 250, 600))))

    gr <- GenomicRanges::makeGRangesFromDataFrame(data.frame(
        seqnames = c("chr1", "chr2"),
        start = c(400, 300),
        end = c(417, 311),
        name = as.character(c("Rank_1", "Rank_2")),
        score = as.integer(c(200, 300)),
        signalValue = as.integer(c(10, 3)),
        pValue = as.integer(c(20, 20)),
        qValue = as.integer(c(20, 19)),
        peak = as.integer(c(930, 220)), stringsAsFactors = FALSE),
        keep.extra.columns = TRUE)

    expect_identical(loadSampleBed("../test_data/with_metadata.bed",
                                   metadata_cols = c("name", "score","strand",
                                                     "signalValue", "pValue",
                                                     "qValue", "peak")), gr)

    expect_warning(loadSampleBed("../test_data/with_metadata.bed"))
    expect_warning(loadSampleBed("../test_data/with_metadata.bed"),
                   metadata_cols = c("name", "score","strand"))
    expect_error(loadSampleBed("tests/test_data/twocol.bed"))

})

test_that("loadBed reads BED files for a mark for the cohort", {

    metadata <- data.frame(Sample = c("A", "B"),
                           TestMark = c("../test_data/basic.bed",
                                        "../test_data/basic.bed"),
                           stringsAsFactors = FALSE)

    exp <- list(A = GRanges(seqnames = rep("chr1", 3),
                            ranges = IRanges(start = c(100, 150, 500),
                                             end = c(200, 250, 600))),
                B = GRanges(seqnames = rep("chr1", 3),
                            ranges = IRanges(start = c(100, 150, 500),
                                             end = c(200, 250, 600))))

    expect_equal(loadBed(metadata, "TestMark"), exp)

    metadata2 <- data.frame(Sample = c("A", "B"),
                           TestMark = c("../test_data/with_metadata.bed",
                                        "../test_data/with_metadata.bed"),
                           stringsAsFactors = FALSE)

    gr <- GenomicRanges::makeGRangesFromDataFrame(data.frame(
        seqnames = c("chr1", "chr2"),
        start = c(400, 300),
        end = c(417, 311),
        name = as.character(c("Rank_1", "Rank_2")),
        score = as.integer(c(200, 300)),
        signalValue = as.integer(c(10, 3)),
        pValue = as.integer(c(20, 20)),
        qValue = as.integer(c(20, 19)),
        peak = as.integer(c(930, 220)), stringsAsFactors = FALSE),
        keep.extra.columns = TRUE)

    expect_equal(loadBed(metadata2, "TestMark",
                         metadata_cols = c("name", "score","strand",
                                           "signalValue", "pValue",
                                           "qValue", "peak")),
                 list(A = gr,
                      B = gr))

    metadata3 <- data.frame(Sample = c("A", "B"),
                            TestMark = c("../test_data/basic.bed",
                                         "../test_data/with_metadata.bed"),
                            stringsAsFactors = FALSE)

    expect_warning(loadBed(metadata3, "TestMark"))

    expect_error(loadBed(metadata3, "Nonexistent"), "mark was not found")

})
