context("Summary strategy for feature matrix construction")

test_that("summarizePeaks generates summary-stats feature matrices", {

    region <- GRanges(seqnames = "chr1",
                      ranges = IRanges(start = 100, end = 600))

    gr <- GenomicRanges::makeGRangesFromDataFrame(data.frame(
        seqnames = c("chr1", "chr1"),
        start = c(150, 400),
        end = c(400, 600),
        name = as.character(c("Rank_1", "Rank_2")),
        score = as.integer(c(200, 300)),
        signalValue = as.integer(c(10, 3)),
        pValue = as.integer(c(20, 20)),
        qValue = as.integer(c(20, 19)),
        peak = as.integer(c(930, 220)), stringsAsFactors = FALSE),
        keep.extra.columns = TRUE)

    lpk <- localPeaks(region = region,
                peaks = list(A = gr, B = gr),
                samples = c("A", "B"))

    summary <- data.frame("score_mean" = c(250, 250),
                            "signalValue_mean" = c(6.5, 6.5),
                            "score_median" = c(250, 250),
                            "signalValue_median" = c(6.5, 6.5),
                            "score_max" = c(300, 300),
                            "signalValue_max" = c(10, 10),
                            "fraction_region_in_peaks" = c(0.9021956,
                                                           0.9021956))
    rownames(summary) <- c("A", "B")
    names(summary) <- paste0("H3K4me3_", names(summary))

    expect_equal(summarizePeaks(lpk, "H3K4me3",
                                cols = c("score", "signalValue")),
                 as.matrix(summary))

    expect_error(summarizePeaks(lpk, "H3K4me3"), "missing")
    expect_error(summarizePeaks(lpk, "H3K4me3", cols = "name"), "not numeric")
    expect_error(summarizePeaks(lpk, "H3K4me3", cols = c("name", "qValue")),
                                "not numeric")

})


test_that("Whole-region strategy accomodates regions with no peaks", {

    empty_pks <- GRanges(seqnames = c(), ranges = IRanges(start = c(), end = c()))
    region <- GRanges(seqnames = "chr2",
                        ranges = IRanges(start = 100, end = 300))

    expect_equal(suppressWarnings(peakOverlap(region, empty_pks)), 0)

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
    mcols(pks$A)$signalValue <- c(5, 10, 15)
    mcols(pks$B)$signalValue <- c(8, 3, 10)

    nopk <- suppressWarnings(retrievePeaks(pks, metadata, region))
    nopk_summary <- data.frame("signalValue_mean" = c(0, 0),
                               "signalValue_median" = c(0, 0),
                               "signalValue_max" = c(0, 0),
                               "fraction_region_in_peaks" = c(0, 0))

    rownames(nopk_summary) <- c("A", "B")
    names(nopk_summary) <- paste0("H3K4me3_", names(nopk_summary))

    expect_equal(summarizePeaks(nopk, "H3K4me3", cols = "signalValue"),
                as.matrix(nopk_summary))

    nopk_summary2 <- nopk_summary %>%
        dplyr::select_("-H3K4me3_fraction_region_in_peaks")

    expect_equal(summarizePeaks(nopk, "H3K4me3",
                                cols = "signalValue",
                                fraction = FALSE),
                 as.matrix(nopk_summary2))

})
