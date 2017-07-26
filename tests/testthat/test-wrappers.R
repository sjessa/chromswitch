context("Chromswitch wrapper functions")


test_that("The whole-region strategy wrapper properly executes the analysis", {

    samples <- c("brain1", "brain2", "brain3", "other1", "other2", "other3")
    outfiles <- system.file("extdata", paste0(samples, ".H3K4me3.bed"),
    package = "chromswitch")
    groups <- c(rep("Brain", 3), rep("Other", 3))

    metadata <- data.frame(Sample = samples,
        H3K4me3 = outfiles,
        Group = groups,
        stringsAsFactors = FALSE)

    regions <- GenomicRanges::GRanges(seqnames = c("chr19", "chr19"),
        ranges = IRanges::IRanges(start = c(54924104, 54874318),
                                  end = c(54929104, 54877536)))

    output <- data.frame(
        region = c("chr19:54924104-54929104", "chr19:54874318-54877536"),
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
        brain1 = as.integer(c(1, 1)),
        brain2 = as.integer(c(1, 1)),
        brain3 = as.integer(c(1, 2)),
        other1 = as.integer(c(2, 1)),
        other2 = as.integer(c(2, 1)),
        other3 = as.integer(c(2, 2)), stringsAsFactors = FALSE
    )

    call <- pryr::partial(callWholeRegion,
                  query = regions,
                  peaks = H3K4me3,
                  metadata = metadata)

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
        brain1 = as.integer(c(1, 1)),
        brain2 = as.integer(c(1, 2)),
        brain3 = as.integer(c(1, 1)),
        other1 = as.integer(c(2, 2)),
        other2 = as.integer(c(2, 2)),
        other3 = as.integer(c(2, 1)), stringsAsFactors = FALSE
    )

    expect_equal(call(normalize = FALSE,
                      mark = "H3K4me3",
                      summarize_columns = c("pValue", "qValue", "signalValue"),
                      heatmap = FALSE),
                 output_nonorm, tolerance = 1e-2)

    expect_error(call(normalize = FALSE,
                      mark = "H3K4me3",
                      summarize_columns = c("pValue", "qValue", "signalValue"),
                      titles = "test"), "one title per query")

})
