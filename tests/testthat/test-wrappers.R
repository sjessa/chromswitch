context("Chromswitch wrapper functions")


test_that("The summary strategy wrapper properly executes the analysis", {
    skip_on_os("windows")

    samples <- c("E068", "E071", "E074", "E101", "E102", "E110")
    bedfiles <- system.file("extdata", paste0(samples, ".H3K4me3.bed"),
                            package = "chromswitch")
    groups <- c(rep("Brain", 3), rep("Other", 3))

    metadata <- data.frame(Sample = samples,
                           H3K4me3 = bedfiles,
                           Condition = groups,
                           stringsAsFactors = FALSE)

    regions <- GenomicRanges::GRanges(seqnames = c("chr19", "chr19"),
                                      ranges = IRanges::IRanges(start = c(54924104, 54874318),
                                                                end = c(54929104, 54877536)))

    mcols(regions)$name <- c("test1", "test2")

    output_0 <- data.frame(
        query = c("chr19:54924104-54929104", "chr19:54874318-54877536"),
        k = as.integer(c(2, 2)),
        Average_Silhouette = c(0.9822687, 0.5608182),
        Consensus = c(1, -0.07207207),
        E068 = as.integer(c(1, 1)),
        E071 = as.integer(c(1, 1)),
        E074 = as.integer(c(1, 2)),
        E101 = as.integer(c(2, 1)),
        E102 = as.integer(c(2, 1)),
        E110 = as.integer(c(2, 2)), stringsAsFactors = FALSE
    )

    output <- data.frame(
        query = c("chr19:54924104-54929104", "chr19:54874318-54877536"),
        k = as.integer(c(2, 2)),
        Average_Silhouette = c(0.911, 0.379),
        Consensus = c(1, 0.156),
        E068 = as.integer(c(1, 1)),
        E071 = as.integer(c(1, 1)),
        E074 = as.integer(c(1, 1)),
        E101 = as.integer(c(2, 1)),
        E102 = as.integer(c(2, 1)),
        E110 = as.integer(c(2, 2)), stringsAsFactors = FALSE
    )

    call <- function(...) {

        callSummary(query = regions,
                    peaks = H3K4me3,
                    metadata = metadata, ...)

    }

    # No summarize_cols, no normalize_cols
    expect_equal(call(mark = "H3K4me3"), output_0, tolerance = 1e-2)

    expect_equal(call(normalize_columns = c("qValue", "pValue", "signalValue"),
                      mark = "H3K4me3",
                      summarize_columns = c("pValue", "qValue", "signalValue"),
                      heatmap = FALSE),
                 output, tolerance = 1e-2)

    # Testing that normalize_cols gets summarize_cols by default
    expect_equal(call(mark = "H3K4me3",
                      summarize_columns = c("pValue", "qValue", "signalValue"),
                      heatmap = FALSE),
                 output, tolerance = 1e-2)

    expect_error(call(mark = "H3K4me3", normalize = FALSE,
                      filter = TRUE,
                      filter_columns = "pValue",
                      filter_thresholds = c(5, 6),
                      summarize_columns = c("pValue", "qValue", "signalValue"),
                      heatmap = FALSE),
                 "one threshold per column")

    expect_error(call(mark = "H3K4me3", normalize = FALSE,
                      filter = TRUE,
                      filter_columns = "pValue",
                      filter_thresholds = NULL,
                      summarize_columns = c("pValue", "qValue", "signalValue"),
                      heatmap = FALSE),
                 "specify thresholds")

    expect_error(call(mark = "H3K4me3", normalize = FALSE,
                      filter = FALSE,
                      summarize_columns = NULL,
                      fraction = FALSE, n = FALSE),
                 "cannot construct")

    expect_error(call(mark = "H3K4me3", normalize = FALSE,
                      filter = TRUE,
                      filter_columns = "pValue",
                      filter_thresholds = NULL,
                      summarize_columns = c("pValue", "qValue", "signalValue"),
                      heatmap = FALSE),
                 "specify thresholds")

    output_nonorm <- data.frame(
        query = c("chr19:54924104-54929104", "chr19:54874318-54877536"),
        k = as.integer(c(2, 2)),
        Average_Silhouette = c(0.695, 0.461),
        Consensus = c(1, 0.0174),
        E068 = as.integer(c(1, 1)),
        E071 = as.integer(c(1, 2)),
        E074 = as.integer(c(1, 1)),
        E101 = as.integer(c(2, 2)),
        E102 = as.integer(c(2, 2)),
        E110 = as.integer(c(2, 1)), stringsAsFactors = FALSE
    )

    expect_equal(call(normalize = FALSE,
                      mark = "H3K4me3",
                      summarize_columns = c("pValue", "qValue", "signalValue"),
                      heatmap = FALSE),
                 output_nonorm, tolerance = 1e-2)

    expect_equal(call(normalize = FALSE,
                      mark = "H3K4me3",
                      summarize_columns = c("pValue", "qValue", "signalValue"),
                      heatmap = TRUE),
                 output_nonorm, tolerance = 1e-2)

    file.remove(paste0(GRangesToCoord(regions[1]), ".pdf"))
    file.remove(paste0(GRangesToCoord(regions[2]), ".pdf"))

    expect_equal(call(normalize = FALSE,
                      mark = "H3K4me3",
                      summarize_columns = c("pValue", "qValue", "signalValue"),
                      heatmap = TRUE,
                      outdir = "outdir"),
                 output_nonorm, tolerance = 1e-2)

    unlink("outdir", recursive = TRUE)

    expect_error(call(normalize = FALSE,
                      mark = "H3K4me3",
                      summarize_columns = c("pValue", "qValue", "signalValue"),
                      titles = "test"), "one title per query")

    output_state <- data.frame(
        query = c("chr19:54924104-54929104", "chr19:54874318-54877536"),
        k = as.integer(c(2, 2)),
        Average_Silhouette = c(0.911, 0.3793),
        Consensus = c(1, 0.1560355),
        state = c("ON", NA),
        E068 = as.integer(c(1, 1)),
        E071 = as.integer(c(1, 1)),
        E074 = as.integer(c(1, 1)),
        E101 = as.integer(c(2, 1)),
        E102 = as.integer(c(2, 1)),
        E110 = as.integer(c(2, 2)), stringsAsFactors = FALSE
    )

    expect_equal(call(normalize_columns = c("qValue", "pValue", "signalValue"),
                      mark = "H3K4me3",
                      summarize_columns = c("pValue", "qValue", "signalValue"),
                      heatmap = FALSE,
                      estimate_state = TRUE,
                      signal_col = "signalValue",
                      test_condition = "Brain"),
                 output_state, tolerance = 1e-2)

    expect_error(call(normalize_columns = c("qValue", "pValue", "signalValue"),
                      mark = "H3K4me3",
                      summarize_columns = c("pValue", "qValue", "signalValue"),
                      heatmap = FALSE,
                      estimate_state = TRUE), "condition")

})


test_that("The binary strategy classifies regions correctly", {
    skip_on_os("windows")

    samples <- c("E068", "E071", "E074", "E101", "E102", "E110")
    bedfiles <- system.file("extdata", paste0(samples, ".H3K4me3.bed"),
                            package = "chromswitch")
    groups <- c(rep("Brain", 3), rep("Other", 3))

    metadata <- data.frame(Sample = samples,
                           H3K4me3 = bedfiles,
                           Condition = groups,
                           stringsAsFactors = FALSE)

    regions <- GenomicRanges::GRanges(seqnames = c("chr19", "chr19"),
                                      ranges = IRanges::IRanges(start = c(54924104, 54892830),
                                                                end = c(54929104, 54897288)))

    output <- data.frame(query = c("chr19:54924104-54929104",
                                   "chr19:54892830-54897288"),
                         k = c(2, 2),
                         Average_Silhouette = c(1, 0),
                         Consensus = c(1.0000000, 0.1560355),
                         E068 = c(1, 1),
                         E071 = c(1, 1),
                         E074 = c(1, 1),
                         E101 = c(2, 1),
                         E102 = c(2, 1),
                         E110 = c(2, 2), stringsAsFactors = FALSE)

    expect_equal(suppressWarnings(callBinary(query = regions,
                                             peaks = H3K4me3,
                                             metadata = metadata,
                                             filter = FALSE,
                                             reduce = TRUE)),
                 output, tolerance = 1e-4)

    output2 <- data.frame(query = c("chr19:54924104-54929104",
                                    "chr19:54892830-54897288"),
                          k = c(2, 2),
                          Average_Silhouette = c(0.8333333, 0),
                          Consensus = c(0.1560355, 0.1560355),
                          E068 = c(1, 1),
                          E071 = c(2, 1),
                          E074 = c(1, 1),
                          E101 = c(1, 1),
                          E102 = c(1, 1),
                          E110 = c(1, 2), stringsAsFactors = FALSE)

    expect_equal(suppressWarnings(callBinary(query = regions,
                                             peaks = H3K4me3,
                                             metadata = metadata,
                                             filter = TRUE,
                                             filter_columns = c("signalValue"),
                                             # A very extreme threshold, to make
                                             # sure the filtering works
                                             filter_thresholds = c(25))),
                 output2, tolerance = 1e-5)

    expect_error(callBinary(query = regions,
                            peaks = H3K4me3,
                            metadata = metadata,
                            filter = TRUE),
                 "provide names of columns to filter")

    output3 <- data.frame(query = "chr19:54924104-54929104",
                          k = 2,
                          Average_Silhouette = 0.4065566,
                          Consensus = 0.1560355,
                          n_features = 2,
                          E068 = 1,
                          E071 = 1,
                          E074 = 2,
                          E101 = 1,
                          E102 = 1,
                          E110 = 1, stringsAsFactors = FALSE)

    expect_equal(callBinary(query = regions[1],
                            peaks = H3K4me3,
                            metadata = metadata,
                            filter = FALSE,
                            heatmap = TRUE,
                            p = 0.9,
                            optimal_clusters = FALSE,
                            n_features = TRUE), output3, tolerance = 1e-5)

    file.remove(paste0(GRangesToCoord(regions[1]), ".pdf"))

})
