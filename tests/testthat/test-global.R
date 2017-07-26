context("Global processing functions")

test_that("winsorNorm normalizes vectors", {
    expect_equal(winsorNorm(seq(0, 10), trim = 0.1),
                 c(0, 0, 0.125, 0.250, 0.375, 0.500, 0.625, 0.75, 0.875, 1, 1))

})

test_that("normalizePeaks normalizes peak stats genome-wide within samples", {

    exp <- list(A = GRanges(seqnames = rep("chr1", 3),
                            ranges = IRanges(start = c(100, 150, 500),
                                             end = c(200, 250, 600))),
                B = GRanges(seqnames = rep("chr1", 3),
                            ranges = IRanges(start = c(100, 150, 500),
                                             end = c(200, 250, 600))))

    expect_error(normalizePeaks(exp))
    expect_error(normalizePeaks(exp, "Nonexistent"),
                 "not found in the metadata")

    gr <- makeGRangesFromDataFrame(data.frame(
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

    cohort <- list(A = gr,
                   B = gr)

    gr_norm <- makeGRangesFromDataFrame(data.frame(
        seqnames = c("chr1", "chr2"),
        start = c(400, 300),
        end = c(417, 311),
        name = as.character(c("Rank_1", "Rank_2")),
        score = as.integer(c(200, 300)),
        signalValue = as.integer(c(1, 0)),
        pValue = as.integer(c(NA, NA)), # Because the vector is c(20, 20)
        qValue = as.integer(c(20, 19)),
        peak = as.integer(c(930, 220)), stringsAsFactors = FALSE),
        keep.extra.columns = TRUE)

    expect_equal(normalizePeaks(cohort, columns = c("signalValue", "pValue")),
                                list(A = gr_norm,
                                     B = gr_norm))

    expect_error(normalizePeaks(cohort, columns = "name"),
                 "non-numeric argument")

    expect_error(normalizePeaks(cohort, columns = "name", tail = 5),
                 "provide a value")

})


test_that("filterPeaks filters peaks according to specified thresholds", {

    gr <- makeGRangesFromDataFrame(data.frame(
        seqnames = c("chr1", "chr2"),
        start = c(400, 300),
        end = c(417, 311),
        name = as.character(c("Rank_1", "Rank_2")),
        score = as.integer(c(200, 300)),
        signalValue = as.integer(c(10, 3)),
        pValue = as.integer(c(20, 20)),
        qValue = as.integer(c(20, 15)),
        peak = as.integer(c(930, 220)), stringsAsFactors = FALSE),
        keep.extra.columns = TRUE)

    gr_filt <- makeGRangesFromDataFrame(data.frame(
        seqnames = "chr1",
        start = 400,
        end = 417,
        name = as.character("Rank_1"),
        score = as.integer(200),
        signalValue = as.integer(10),
        pValue = as.integer(20),
        qValue = as.integer(20),
        peak = as.integer(930), stringsAsFactors = FALSE),
        keep.extra.columns = TRUE)

    GenomeInfoDb::seqlevels(gr_filt) <- c("chr1", "chr2")

    expect_identical(filterPeaks(list(A = gr, B = gr),
                                columns = c("signalValue", "qValue"),
                                thresholds = c(5, 17)),
                 list(A = gr_filt, B = gr_filt))

    expect_error(filterPeaks(list(A = gr, B = gr),
                                columns = c("signalValue", "qValue"),
                                thresholds = 5), "threshold per column")

    expect_error(filterPeaks(list(A = gr, B = gr),
                            columns = c("signalValue", "qValue"),
                            thresholds = c(5, 6, 7)), "threshold per column")

    expect_error(filterPeaks(list(A = gr, B = gr),
                            columns = c("name", "qValue"),
                            thresholds = c(5, 6)), "not numeric")

    expect_error(filterPeaks(list(A = gr, B = gr),
                            columns = c("nonexistent", "qValue"),
                            thresholds = c(5, 6)), "not found")

})

