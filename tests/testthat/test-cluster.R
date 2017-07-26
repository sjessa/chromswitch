context("Hierarchical clustering over the feature-by-sample matrix")


test_that("Hierarchical clustering finds clusters from feature matrix", {

    samples <- c("brain1", "brain2", "brain3", "other1", "other2", "other3")
    outfiles <- system.file("extdata", paste0(samples, ".H3K4me3.bed"),
    package = "chromswitch")
    groups <- c(rep("Brain", 3), rep("Other", 3))

    metadata <- data.frame(Sample = samples,
        H3K4me3 = outfiles,
        Group = groups,
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
        Consensus_top = 1,
        brain1 = 1, brain2 = 1, brain3 = 1,
        other1 = 2, other2 = 2, other3 = 2, stringsAsFactors = FALSE)

    expect_equal(cluster(ft_mat, metadata, region, heatmap = FALSE),
                cluster_out)

    expect_equal(cluster(ft_mat, metadata, region),
                 cluster_out)

    # Clean up
    file.remove(paste0(GRangesToCoord(region), ".pdf"))

})
