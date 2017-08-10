context("Hierarchical clustering over the feature-by-sample matrix")

test_that("The ranking method correctly estimates state/trajectory", {


    a <- data.frame(Condition = c("C1", "C1", "C2"),
                    fc = c(10, 9, 8))

    b <- data.frame(Condition = c("C2", "C1", "C2"),
                    fc = c(10, 9, 8))

    c <- data.frame(Condition = c("C1", "C2", "C1"),
                    fc = c(10, 9, 8))

    d <- data.frame(Condition = c("C2", "C1", "C1"),
                    fc = c(10, 9, 8))

    expect_equal(estimateState(a), "ON")
    expect_equal(estimateState(b), NA)
    expect_equal(estimateState(c), NA)
    expect_equal(estimateState(d), "OFF")

})


test_that("The clustering function returns the correct set of clusters", {

    ft_mat <- data.frame("ft1" = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE),
                         "ft2" = c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE))

    expect_equal(clusterValidityPerK(ft_mat),
                 data.frame(k = c(2, 3, 4, 5),
                            Average_Silhouette = c(0.6206429, 1.0, 2/3, 1/3)),
                 tolerance = 1e-3)

    best_k <- dplyr::tibble(k = as.integer(3), Average_Silhouette = 1)
    expect_equal(getK(ft_mat, optimal_clusters = TRUE), best_k)

    samples <- c("E068", "E071", "E074", "E101", "E102", "E110")
    outfiles <- system.file("extdata", paste0(samples, ".H3K4me3.bed"),
                            package = "chromswitch")
    groups <- c(rep("Brain", 3), rep("Other", 3))

    metadata <- data.frame(Sample = samples,
                           H3K4me3 = outfiles,
                           Condition = groups,
                           stringsAsFactors = FALSE)

    ft_mat <- data.frame(
        H3K4me3_qValue_mean = c(46, 66, 20, 22, 0, 0),
        H3K4me3_signalValue_mean = c(15, 17, 5, 4, 0, 0),
        H3K4me3_qValue_max = c(93, 115, 47, 50, 0, 0),
        H3K4me3_signalValue_max = c(22, 26, 8, 7, 0, 0),
        H3K4me3_fraction_region_in_peaks = c(0.5, 0.5, 0.2, 0.21, 0, 0)
    )

    rownames(ft_mat) <- metadata$Sample

    region = GenomicRanges::GRanges(seqnames = "chr19",
                                    ranges = IRanges::IRanges(start = 54924104,
                                                              end = 54929104))

    cluster_out <- data.frame(
        region = GRangesToCoord(region),
        k = 3,
        Average_Silhouette = 0.8231449,
        Purity = 0.8333333,
        Entropy = 0.5793802,
        ARI = 0.2424242,
        NMI = 0.5295406,
        Homogeneity = 0.6666667,
        Completeness = 0.4206198,
        V_measure = 0.5158037,
        Consensus = 0.4292562,
        E068 = 1, E071 = 1, E074 = 2,
        E101 = 2, E102 = 3, E110 = 3, stringsAsFactors = FALSE)

    expect_equal(cluster(ft_mat, metadata, region, optimal_clusters = TRUE),
                 cluster_out, tolerance = 1e-6)

    cluster_out_not_optimal <- data.frame(
        region = GRangesToCoord(region),
        k = 2,
        Average_Silhouette = 0.6145993,
        Purity = 0.8333333,
        Entropy = 0.4591479,
        ARI = 0.3243243,
        NMI = 0.4791388,
        Homogeneity = 0.4591479,
        Completeness =  0.5,
        V_measure = 0.478704,
        Consensus = 0.427389,
        E068 = 1, E071 = 1, E074 = 2,
        E101 = 2, E102 = 2, E110 = 2, stringsAsFactors = FALSE)

    expect_equal(cluster(ft_mat, metadata, region, optimal_clusters = FALSE),
                 cluster_out_not_optimal, tolerance = 1e-6)

    ft_mat <- data.frame(no_peak = c(rep(TRUE, length(metadata$Sample))))
    rownames(ft_mat) <- metadata$Sample

    cluster_out2 <- data.frame(
        region = GRangesToCoord(region),
        k = 2,
        Average_Silhouette = 0,
        Purity = 0.6666667,
        Entropy = 0.4591479,
        ARI = 0,
        NMI = 0.2367466,
        Homogeneity = 0.1908745,
        Completeness = 0.293643,
        V_measure = 0.2313599,
        Consensus = 0.1560355,
        E068 = 1, E071 = 1, E074 = 1,
        E101 = 1, E102 = 1, E110 = 2, stringsAsFactors = FALSE)
})


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
        H3K4me3_qValue_mean = c(46, 66, 27, 0, 0, 0),
        H3K4me3_signalValue_mean = c(15, 17, 11, 0, 0, 0),
        H3K4me3_qValue_max = c(73, 115, 47, 0, 0, 0),
        H3K4me3_signalValue_max = c(22, 26, 17, 0, 0, 0),
        H3K4me3_fraction_region_in_peaks = c(0.5, 0.5, 0.5, 0, 0, 0)
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

    cluster_out3 <- data.frame(
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
        state = "ON",
        n_features = 5,
        E068 = 1, E071 = 1, E074 = 1,
        E101 = 2, E102 = 2, E110 = 2, stringsAsFactors = FALSE)

    expect_equal(cluster(ft_mat, metadata, region,
                         n_features = TRUE,
                         estimate_state = TRUE,
                         signal_col = "signalValue",
                         mark = "H3K4me3",
                         test_condition = "Brain"),
                 cluster_out3)

    expect_error(cluster(ft_mat, metadata, region,
                         n_features = TRUE,
                         estimate_state = TRUE), "signal value")


    expect_error(cluster(ft_mat, metadata, region,
                         n_features = TRUE,
                         estimate_state = TRUE,
                         signal_col = "H3K4me3"), "condition")

    expect_error(cluster(ft_mat, metadata, region,
                         n_features = TRUE,
                         estimate_state = TRUE,
                         signal_col = "H3K4me3",
                         test_condition = "Brain"), "mark")

})

test_that("Clustering function assigns number of features", {

    region2 <- GRanges(seqnames = "chr1",
                       ranges = IRanges(start = 100, end = 300))

    pks2 <- GRanges(seqnames = c(), ranges = IRanges(start = c(), end = c()))
    lp2 <- localPeaks(region2, list(A = pks2, B = pks2, C = pks2, D = pks2),
                      c("A", "B", "C", "D"))

    meta <- data.frame(Sample = c("A", "B", "C", "D"),
                       Condition = c("Brain", "Brain", "Other", "Other"))

    ft_mat <- suppressWarnings(binarizePeaks(lp2, 0.5))

    expect_equal(dplyr::select(suppressWarnings(cluster(ft_mat, meta, region2,
                                                 n_features = TRUE)),
                        "n_features"),
                 data.frame(n_features = 0, stringsAsFactors = FALSE))

})
