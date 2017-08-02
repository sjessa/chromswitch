context("Hierarchical clustering over the feature-by-sample matrix")


test_that("Hierarchical clustering finds clusters from feature matrix", {

    samples <- c("E068", "E071", "E074", "E101", "E102", "E110")
    outfiles <- system.file("extdata", paste0(samples, ".H3K4me3.bed"),
    package = "chromswitch")
    groups <- c(rep("Brain", 3), rep("Other", 3))

    metadata <- data.frame(Sample = samples,
        H3K4me3 = outfiles,
        Condition = groups,
        stringsAsFactors = FALSE)

    ft_mat <- data.frame(
        qValue_mean = c(46, 66, 27, 0, 0, 0),
        signalValue_mean = c(15, 17, 11, 0, 0, 0),
        qValue_max = c(73, 115, 47, 0, 0, 0),
        signalValue_max = c(22, 26, 17, 0, 0, 0),
        fraction_region_in_peaks = c(0.5, 0.5, 0.5, 0, 0, 0)
    )

    rownames(ft_mat) <- metadata$Sample

    region = GenomicRanges::GRanges(seqnames = "chr19",
                                    ranges = IRanges::IRanges(start = 54924104,
                                                            end = 54929104))

    cluster_out <- data.frame(
        region = GRangesToCoord(region),
        k = 2,
        Average_Silhouette = 0.6883056,
        Purity = 1,
        Entropy = 0,
        ARI = 1,
        NMI = 1,
        Homogeneity = 1,
        Completeness = 1,
        V_measure = 1,
        Consensus = 1,
        E068 = 1, E071 = 1, E074 = 1,
        E101 = 2, E102 = 2, E110 = 2, stringsAsFactors = FALSE)

    expect_equal(cluster(ft_mat, metadata, region),
                cluster_out)

    expect_equal(cluster(ft_mat, metadata, region, heatmap = TRUE),
                 cluster_out)

    expect_equal(cluster(ft_mat, metadata, region, heatmap = TRUE,
                        outdir = "."),
                 cluster_out)

    # Clean up
    file.remove(paste0(GRangesToCoord(region), ".pdf"))

    region2 <- GenomicRanges::GRanges(seqnames = "chr19",
                                    ranges = IRanges::IRanges(start = 54924104,
                                                              end = 54929104))

    mcols(region2)$name <- "Test"

    cluster_out2 <- data.frame(
        region = GRangesToCoord(region),
        name = "Test",
        k = 2,
        Average_Silhouette = 0.6883056,
        Purity = 1,
        Entropy = 0,
        ARI = 1,
        NMI = 1,
        Homogeneity = 1,
        Completeness = 1,
        V_measure = 1,
        Consensus = 1,
        E068 = 1, E071 = 1, E074 = 1,
        E101 = 2, E102 = 2, E110 = 2, stringsAsFactors = FALSE)

    expect_equal(cluster(ft_mat, metadata, region2, heatmap = FALSE),
                 cluster_out2)

})
