context("Local processing functions")

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
