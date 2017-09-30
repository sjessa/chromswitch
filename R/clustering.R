# ---------------------------------------------------------------------------- #
#
# Clustering function: takes in a feature matrix and returns clusters
# scores in order to call a chromatin switch
#
# ---------------------------------------------------------------------------- #


#' cluster
#'
#' Given a sample-by-feature matrix and sample-associated metadata including
#' their biological condition groupings, cluster samples hierarchically
#' and use external cluster validity measures (Adjusted Rand Index,
#' Normalized Mutual Information, and V measure) to assess the agreement between
#' the inferred clusters and the biological conditions. Optionally, produce
#' a heatmap reflecting the hierarchical clustering result.
#'
#' @param ft_mat matrix where columns are features and rows are samples as
#' returned by \code{\link{summarizePeaks}} or \code{\link{binarizePeaks}}
#' @param metadata A dataframe with a column "Sample" which stores
#' the sample identifiers, and at least one column, "Condition", which stores
#' the biological condition labels of the samples
#' @param region GRanges object specifying the query region
#' @param heatmap (Optional) Logical value indicating whether to plot
#' the heatmap for hierarchical clustering. Default: FALSE
#' @param title (Optional) If \code{heatmap} is TRUE, specify the title of the
#' plot, which will also be used for the output file name in PDF format
#' @param outdir (Optional) String specifying the name of the directory where
#' PDF of heatmaps should be saved
#' @param optimal_clusters (Optional) Logical value indicate whether to cluster
#' samples into two groups, or to find the optimal clustering solution by
#' choosing the set of clusters which maximizes the Average Silhouette width
#' @param n_features (Optional) Logical value indicating whether to include
#' a column "n_features" in the output storing the number of features in the
#' feature matrix constructed for the region, which may be useful for
#' understanding the behaviour of the binary strategy for constructing
#' feature matrices. Default: FALSE
#' @param estimate_state (Optional) Logical value indicating whether to include
#' a column "state" in the output specifying the estimated chromatin state of
#' a test condition. The state will be on of "ON", "OFF", or NA, where the
#' latter results if a binary switch between the conditions is unclear.
#' Default: FALSE.
#' @param method (Optional) If \code{estimate_state} is TRUE, one of "summary"
#' or "binary", specifying which method was used to construct the feature
#' matrix in \code{ft_mat}
#' @param test_condition (Optional) If \code{estimate_state} is TRUE, string
#' specifying one of the two biological condtions in \code{metadata$Condition}
#' for which to estimate chromatin state.
#' @param signal_col (Optional) If \code{estimate_state} is TRUE, and
#' \code{method} is "summary", string
#' specifying the name of the column in the original peak files which
#' corresponds to the level of enrichment in the region, e.g. fold change
#' @param mark (Optional) If \code{estimate_state} is TRUE, and \code{method}
#' is "summary",string specifying
#' the name of the mark for which \code{ft_mat} was constructed
#'
#' @examples
#' samples <- c("E068", "E071", "E074", "E101", "E102", "E110")
#' outfiles <- system.file("extdata", paste0(samples, ".H3K4me3.bed"),
#' package = "chromswitch")
#' Conditions <- c(rep("Brain", 3), rep("Other", 3))
#'
#' metadata <- data.frame(Sample = samples,
#'     H3K4me3 = outfiles,
#'     Condition = Conditions,
#'     stringsAsFactors = FALSE)
#'
#' region <- GRanges(seqnames = "chr19",
#'     ranges = IRanges(start = 54924104, end = 54929104))
#'
#' lpk <- retrievePeaks(H3K4me3,
#'     metadata = metadata,
#'     region = region)
#'
#' ft_mat <- summarizePeaks(lpk, mark = "H3K4me3",
#' cols = c("qValue", "signalValue"))
#'
#' cluster(ft_mat, metadata, region)
#'
#' # Estimate the state of the test condition, "Brain"
#' cluster(ft_mat, metadata, region,
#'     estimate_state = TRUE,
#'     method = "summary",
#'     signal_col = "signalValue",
#'     mark = "H3K4me3",
#'     test_condition = "Brain")
#'
#' @return A dataframe with the region, the number of clusters inferred,
#' the cluster validity statistics, and the cluster assignments for each sample
#'
#' @export
cluster <- function(ft_mat, metadata, region,
                    heatmap = FALSE, title = NULL, outdir = NULL,
                    optimal_clusters = FALSE,
                    n_features = FALSE,
                    estimate_state = FALSE,
                    method = NULL,
                    test_condition = NULL,
                    signal_col = NULL,
                    mark = NULL) {

    if (is(region, "GRangesList")) region <- unlist(region)

    # If only one feature, can't draw a heatmap
    if (ncol(ft_mat) == 1) heatmap = FALSE
    else heatmap = heatmap

    features <- attr(ft_mat, "features")
    ft_mat <- data.matrix(ft_mat)

    if (isTRUE(heatmap)) {

        palette <- grDevices::colorRampPalette(
            c("dodgerblue4", "white", "red"))(n = 100)
        hclust.fun <- function(i) hclust(i, method = "complete")

        conditions <- unique(metadata$Condition)
        conditions_colours <- metadata$Condition

        conditions_colours[conditions_colours == conditions[1]] <-"mediumorchid"
        conditions_colours[conditions_colours == conditions[2]] <- "limegreen"

        if (is.null(title)) title <- GRangesToCoord(region)

        outfile <- ifelse(!is.null(outdir),
                        paste0(outdir, "/", title, ".pdf"),
                        paste0(title, ".pdf"))

        grDevices::pdf(outfile)
        results <- gplots::heatmap.2(ft_mat,
                            dendrogram = "row",
                            trace = "none",
                            col = palette,
                            hclustfun = hclust.fun,
                            RowSideColors = conditions_colours,
                            ylab = "Sample",
                            xlab = "Feature",
                            labCol = colnames(ft_mat),
                            labRow = rownames(ft_mat),
                            key = TRUE,
                            cexCol = 0.8,
                            cexRow = 1.0,
                            margins = c(12, 5),
                            main = title)

        graphics::legend("topright",
                title = "Condition",
                legend = conditions,
                fill = c("mediumorchid", "limegreen"),
                border = c("mediumorchid", "limegreen"),
                cex = 0.8,
                box.lwd = 0,
                bty = "n")
        grDevices::dev.off()

    }

    # Choose the clustering partition with the highest average Silhouette width
    stats <- getK(ft_mat, optimal_clusters = optimal_clusters)

    # Get the clusters
    d_mat <- dist(ft_mat)
    tree <- hclust(d_mat)
    clusters <- cutree(tree, k = stats$k)

    # Calculate external cluster validity stats for the partition by comparing
    # the clusters to the condition labels
    contingency        <- contingencyTable(clusters, metadata)
    stats$Purity       <- purity(contingency)
    stats$Entropy      <- NMF::entropy(contingency)
    stats$ARI          <- mclust::adjustedRandIndex(clusters,
                                                    metadata$Condition)
    stats$NMI          <- NMI(clusters, metadata$Condition)
    stats$Homogeneity  <- homogeneity(contingency)
    stats$Completeness <- completeness(contingency)
    stats$V_measure    <- vMeasure(contingency)
    stats$Consensus    <- mean(x = c(stats$ARI, stats$NMI, stats$V_measure))

    stats <- stats %>% dplyr::select(-c(Purity, Entropy, ARI, NMI,
                                       Homogeneity, Completeness, V_measure))

    clusters_df <- clusters %>%
        as.list() %>%
        as.data.frame(stringsAsFactors = FALSE)
    names(clusters_df) <- names(clusters)

    coord <- GRangesToCoord(region)

    meta_cols <- mcols(region)

    region_df <- data.frame(region = coord,
                            meta_cols,
                            stringsAsFactors = FALSE)

    # Experimental: assign each cluster to a condition and return
    # mean signal (e.g. fold change) in each condition according to the clusters
    # which allows for an estimation of the "trajectory" of the switch
    if (estimate_state) {

        if (!(method %in% c("summary", "binary")))
            stop("If estimate_trajectory is TRUE, specify the method",
                    "used to construct feature matrices, one of 'summary' or",
                    "'binary'.")

        if (is.null(test_condition))
            stop("If estimate_trajectory = TRUE, please specify the name
                of the condition for which to call chromatin state.")

        # 1. Assign each cluster to one of the two conditions
        clust_id <- clusters_df %>%
            t() %>%
            as.data.frame() %>%
            dplyr::mutate(Sample = rownames(.)) %>%
            dplyr::rename_(Cluster = "V1") %>%
            dplyr::inner_join(metadata, by = "Sample") %>%
            dplyr::select(Sample, Condition, Cluster) %>%
            dplyr::mutate(Condition = ifelse(Condition == test_condition,
                                            "C1", "C2"))

        if (method == "summary") {

            if (is.null(signal_col))
                stop("If estimate_trajectory = TRUE, please specify the name
                of the column corresponding to the signal value.")


            if (is.null(mark))
                stop("If estimate_trajectory = TRUE, please specify the name
                of the mark corresponding to the data in ft_mat.")

            col <- ifelse(signal_col == "fraction",
                            paste0(mark, "_fraction_region_in_peaks"),
                            paste0(mark, "_", signal_col, "_mean"))
            ft_mat2 <- ft_mat


        } else if (method == "binary") {

            ft_mat_df <- as.data.frame(ft_mat)

            toLength <- function(i) {

                length <- width(features[i])
                ft_mat_df[,i] %>%
                    dplyr::recode(`1` = length, `0` = as.integer(0)) %>%
                    data.frame
            }

            ft_mat2 <- lapply(seq_along(ft_mat_df), toLength) %>%
                dplyr::bind_cols()
            ft_mat2$olap <- rowSums(ft_mat2)
            ft_mat2 <- ft_mat2 %>% dplyr::mutate(
                frac = olap / (width(region) + 1))
            rownames(ft_mat2) <- rownames(ft_mat)

            col <- "frac"
        }

        # 2. Get the mean signal or fraction of overlap of each cluster
        clust_ft_mat <- ft_mat2 %>% as.data.frame %>%
            dplyr::select_(col) %>%
            dplyr::mutate(Sample = rownames(.)) %>%
            dplyr::inner_join(clust_id, by = "Sample") %>%
            dplyr::group_by(Cluster, Condition) %>%
            dplyr::summarize_at(col, mean, na.rm = TRUE) %>%
            dplyr::arrange_(paste0("desc(", col, ")"))

        # 3. Guess the state of the test condition in the region
        stats$state <- estimateState(clust_ft_mat)
    }

    if (n_features) stats$n_features <- ifelse("no_peak" %in% colnames(ft_mat),
                                                0,
                                                ncol(ft_mat))

    results <- dplyr::bind_cols(region_df, stats, clusters_df)

    return(as.data.frame(results))

}



#' @keywords internal
isC1AtTop <- function(clust_ft_mat) {

    state = FALSE
    i = 1
    current = clust_ft_mat[i, "Condition"]

    # Keep going down the list while the condition is brain
    while(current == "C1") {
        state = TRUE
        i <- i + 1
        current <- clust_ft_mat[i, "Condition"]
    }

    # When it's not C1 anymore, check if there are any C1 clusters below
    leftover <- clust_ft_mat[i:nrow(clust_ft_mat), "Condition"] %>% unlist()
    if ("C1" %in% leftover) state = FALSE

    return(state)
}


#' @keywords internal
isC1AtBottom <- function(clust_ft_mat) {

    # Ask if C1 is at the top when the clusters are ranked from
    # lowest mean FC to highest
    isC1AtTop(clust_ft_mat[order(nrow(clust_ft_mat):1),])
}


#' @keywords internal
estimateState <- function(clust_ft_mat) {

    if (isC1AtTop(clust_ft_mat)) return("ON")
    else if (isC1AtBottom(clust_ft_mat)) return("OFF")
    else return(NA)
}

