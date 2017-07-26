# ---------------------------------------------------------------------------- #
#
# Clustering function: takes in a feature matrix and returns clusters
# scores in order to call a
#
# ---------------------------------------------------------------------------- #


#' cluster
#'
#' @param ft_mat matrix where columns are features and rows are samples as
#' returned by \code{\link{summarizePeaks}} or \code{\link{binarizePeaks}}
#' @param metadata A dataframe with a column "Sample" which stores
#' the sample identifiers, and at least one column, "Group", which stores
#' the biological condition labels of the samples
#' @param region GRanges object specifying the query region
#' @param heatmap (Optional) Logical value indicating whether or not to plot
#' the heatmap for hierarchical clustering. Default: FALSE
#' @param title (Optional) If \code{heatmap} is TRUE, specify the title of the
#' plot, which will also be used for the output file name in PDF format
#' @param outdir Optional, the name of the directory where heatmaps should
#' be saved
#'
#' @examples
#' samples <- c("brain1", "brain2", "brain3", "other1", "other2", "other3")
#' outfiles <- system.file("extdata", paste0(samples, ".H3K4me3.bed"),
#' package = "chromswitch")
#' groups <- c(rep("Brain", 3), rep("Other", 3))
#'
#' metadata <- data.frame(Sample = samples,
#'     H3K4me3 = outfiles,
#'     Group = groups,
#'     stringsAsFactors = FALSE)
#'
#' region <- GenomicRanges::GRanges(seqnames = "chr19",
#'     ranges = IRanges::IRanges(start = 54924104, end = 54929104))
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
#' @return A dataframe with the region, the number of clusters inferred,
#' the cluster validity statistics, and the cluster assignments for each sample
#'
#' @export
cluster <- function(ft_mat, metadata, region,
                    heatmap = FALSE, title = NULL, outdir = NULL) {

    # If only one feature, can't draw a heatmap
    if (ncol(ft_mat) == 1) heatmap = FALSE

    ft_mat <- data.matrix(ft_mat)

    palette <- grDevices::colorRampPalette(c("dodgerblue4",
                                            "white", "red"))(n = 100)
    hclust.fun <- function(i) hclust(i, method = "complete")

    conditions <- unique(metadata$Group)
    conditions_colours <- metadata$Group

    conditions_colours[conditions_colours == conditions[1]] <- "mediumorchid"
    conditions_colours[conditions_colours == conditions[2]] <- "limegreen"

    if (isTRUE(heatmap)) {

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
    stats <- ft_mat %>%
        clusterValidityPerK() %>%
        dplyr::arrange_("desc(Average_Silhouette)") %>%
        dplyr::slice(1)

    d_mat <- dist(ft_mat)
    tree <- hclust(d_mat)
    clusters <- cutree(tree, k = stats$k)

    # Calculate external cluster validity stats for the partition by comparing
    # the clusters to the condition labels
    contingency        <- contingencyTable(clusters, metadata)
    stats$Purity       <- purity(contingency)
    stats$Entropy      <- NMF::entropy(contingency)
    stats$ARI          <- mclust::adjustedRandIndex(clusters, metadata$Group)
    stats$NMI          <- NMI(clusters, metadata$Group)
    stats$Homogeneity  <- homogeneity(contingency)
    stats$Completeness <- completeness(contingency)
    stats$V_measure    <- vMeasure(contingency)
    stats$Consensus_top <- mean(x = c(stats$ARI, stats$NMI, stats$V_measure))

    clusters_df <- clusters %>%
        as.list() %>%
        as.data.frame(stringsAsFactors = FALSE)
    names(clusters_df) <- names(clusters)

    region_df <- data.frame(region = GRangesToCoord(region),
                            stringsAsFactors = FALSE)

    results <- dplyr::bind_cols(region_df, stats, clusters_df)

    return(as.data.frame(results))

}
