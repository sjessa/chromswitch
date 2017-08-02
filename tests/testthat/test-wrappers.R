context("Chromswitch wrapper functions")


test_that("The whole-region strategy wrapper properly executes the analysis", {

    samples <- c("E068", "E071", "E074", "E101", "E102", "E110")
    outfiles <- system.file("extdata", paste0(samples, ".H3K4me3.bed"),
    package = "chromswitch")
    groups <- c(rep("Brain", 3), rep("Other", 3))

    metadata <- data.frame(Sample = samples,
        H3K4me3 = outfiles,
    Condition = groups,
        stringsAsFactors = FALSE)

    regions <- GenomicRanges::GRanges(seqnames = c("chr19", "chr19"),
        ranges = IRanges::IRanges(start = c(54924104, 54874318),
                                  end = c(54929104, 54877536)))

    mcols(regions)$name <- c("test1", "test2")

    output <- data.frame(
        region = c("chr19:54924104-54929104", "chr19:54874318-54877536"),
        name = c("test1", "test2"),
        k = c(2, 2),
        Average_Silhouette = c(0.911, 0.457),
        Purity = c(1, 0.5),
        Entropy = c(0, 0.918),
        ARI = c(1, -0.216),
        NMI = c(1, 0),
        Homogeneity = c(1, 0),
        Completeness = c(1, 0),
        V_measure = c(1, 0),
        Consensus_top = c(1, -0.0721),
        E068 = as.integer(c(1, 1)),
        E071 = as.integer(c(1, 1)),
        E074 = as.integer(c(1, 2)),
        E101 = as.integer(c(2, 1)),
        E102 = as.integer(c(2, 1)),
        E110 = as.integer(c(2, 2)), stringsAsFactors = FALSE
    )

    call <- function(...) {

        callWholeRegion(query = regions,
                       peaks = H3K4me3,
                       metadata = metadata, ...)

    }

    expect_equal(call(normalize_columns = c("qValue", "pValue", "signalValue"),
                    mark = "H3K4me3",
                    summarize_columns = c("pValue", "qValue", "signalValue"),
                    heatmap = FALSE),
                output, tolerance = 1e-2)

    expect_error(call(mark = "H3K4me3",
                    summarize_columns = c("pValue", "qValue", "signalValue"),
                    heatmap = FALSE),
                "provide names of columns to normalize")

    expect_error(call(mark = "H3K4me3", normalize = FALSE,
                    filter = TRUE,
                    summarize_columns = c("pValue", "qValue", "signalValue"),
                    heatmap = FALSE),
                 "provide names of columns to filter")

    expect_error(call(mark = "H3K4me3", normalize = FALSE,
                      filter = TRUE,
                      filter_columns = "pValue",
                      filter_thresholds = c(5, 6),
                      summarize_columns = c("pValue", "qValue", "signalValue"),
                      heatmap = FALSE),
                 "one threshold per column")

    output_nonorm <- data.frame(
        region = c("chr19:54924104-54929104", "chr19:54874318-54877536"),
        name = c("test1", "test2"),
        k = c(2, 2),
        Average_Silhouette = c(0.6949372, 0.4614938),
        Purity = c(1.0000000, 0.6666667),
        Entropy = c(0.0000000, 0.9182958),
        ARI = c(1, -0.1111111),
        NMI = c(1, 0.08170417),
        Homogeneity = c(1, 0.08170417),
        Completeness = c(1, 0.08170417),
        V_measure = c(1, 0.08170417),
        Consensus_top = c(1, 0.01743241),
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

    expect_error(call(normalize = FALSE,
                      mark = "H3K4me3",
                      summarize_columns = c("pValue", "qValue", "signalValue"),
                      titles = "test"), "one title per query")

})


test_that("The position-aware strategy properly executes the analysis", {

    samples <- c("E068", "E071", "E074", "E101", "E102", "E110")
    outfiles <- system.file("extdata", paste0(samples, ".H3K4me3.bed"),
                            package = "chromswitch")
    groups <- c(rep("Brain", 3), rep("Other", 3))

    metadata <- data.frame(Sample = samples,
                           H3K4me3 = outfiles,
                           Condition = groups,
                           stringsAsFactors = FALSE)

    regions <- GenomicRanges::GRanges(seqnames = c("chr19", "chr19"),
                    ranges = IRanges::IRanges(start = c(54924104, 54892830),
                                                end = c(54929104, 54897288)))

    output <- data.frame(region = c("chr19:54924104-54929104",
                                    "chr19:54892830-54897288"),
                        k = c(2, 2),
                        Average_Silhouette = c(1, 0),
                        Purity = c(1, 0.6666667),
                        Entropy = c(0, 0.4591479),
                        ARI = c(1, 0),
                        NMI = c(1, 0.2367466),
                        Homogeneity = c(1, 0.1908745),
                        Completeness = c(1, 0.293643),
                        V_measure = c(1, 0.2313599),
                        Consensus_top = c(1.0000000, 0.1560355),
                        E068 = c(1, 1),
                        E071 = c(1, 1),
                        E074 = c(1, 1),
                        E101 = c(2, 1),
                        E102 = c(2, 1),
                        E110 = c(2, 2), stringsAsFactors = FALSE)

    expect_equal(suppressWarnings(callPositionAware(query = regions,
                                                    peaks = H3K4me3,
                                                    metadata = metadata,
                                                    filter = FALSE,
                                                    reduce = TRUE)),
                output, tolerance = 1e-4)


})
