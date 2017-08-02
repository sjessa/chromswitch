# ---------------------------------------------------------------------------- #
#
# Wrapper functions:
# These combine all the steps in the analysis in order to provide a simple
# way to apply the method to a genomic region for one mark
#
# ---------------------------------------------------------------------------- #


#' callWholeRegion
#'
#' One of two main functions in the \code{chromswitch} package, this function
#' detects a switch in chromatin state in one or
#' more regions given ChIP-seq peak calls for one mark, executing the entire
#' algorithm from preprocessing to evaluating the clustering results,
#' using the whole-region strategy.
#'
#' This strategy constructs a sample-by-feature matrix to use as input for
#' hierarchical clustering by computing, for each sample, a vector of summary
#' statistics based on that sample's peaks in the query region. The summary
#' statistics are generally based on the enrichment statistics associated with
#' each peak as returned by the peak calling too, which might include, for
#' example, a p value and fold change.
#'
#' @param query GRanges list containing one or more genomic regions of interest
#' in which to call a switch. The output dataframe will contain one row per
#' region in \code{query}.
#' @param peaks List of GRanges objects storing peak calls for each sample,
#' where element names correspond to sample IDs
#' @param metadata A dataframe with at least two columns: "Sample" which stores
#' the sample IDs, "Condition", which stores the biological condition labels
#' of the samples
#' @param mark Character specifying the histone mark or ChIP-target,
#' for example, "H3K4me3"
#' @param filter Optional: logical value, filter peaks based on thresholds on
#' peak statistics? Default: FALSE. The filter step is described in
#' \code{\link{filterPeaks}}.
#' @param filter_columns If \code{filter} is TRUE, a chracter vector
#' corresponding to names of columns in the peak metadata by which to filter
#' peaks. If \code{filter} is FALSE, not used.
#' @param filter_thresholds If \code{filter} is TRUE, a numeric vector
#' corresponding to lower cutoffs applied to metadata columns in order to filter
#' peaks. Provide one per column specified in \code{filter_columns}, in the same
#' order. If \code{filter} is FALSE, not used.
#' @param normalize Optional: logical value, normalize peak statistics
#' genome-wide for each sample? Default: TRUE. The normalization step is
#' described in \code{\link{normalizePeaks}}.
#' @param normalize_columns If \code{normalize} is TRUE, a character vector
#' corresponding to names of columns in the peak metadata to normalize
#' genome-wide for each sample. If \code{normalize} is FALSE, not used.
#' @param tail Optional: if \code{normalize} is TRUE, specifies the fraction
#' of extreme values in each tail to bound during normalization. More details at
#' \code{\link{normalizePeaks}}.
#' @param summarize_columns Character vector of column names on which to compute
#' summary statistics during feature matrix construction. These statistics
#' become the features of the matrix.
#' @param length Optional: Logical value, during feature matrix construction,
#' compute the mean, median, and max of peak length? Default: FALSE
#' @param fraction Optional: Logical value, during feature matrix construction,
#' compute the fraction of the region overlapped by peaks? Default: TRUE
#' @param n Optional: Logical value, during feature matrix construction,
#' compute the number of peaks in the region? Default: FALSE
#' @param heatmap Optional: Logical value, plot the heatmap corresponding to
#' the hierarchical clustering result? Default: TRUE
#' @param titles Optional:  if \code{heatmap} is TRUE, a character vector
#' of the same length as \code{query}, specifying the title to use when plotting
#' each heatmap (e.g. a gene name), also reused as the
#' prefix of the name of the file where the heatmap is saved. By default, the
#' title is the genomic coordinates of the region in the form "chrN:start-end"
#' @param outdir Optional: if \code{heatmap} is TRUE, the name of the directory
#' where heatmaps should be saved
#'
#' @return Data frame with one row per region in \code{query}. Contains the
#' coordinates of the region, the number of inferred clusters, the computed
#' cluster validity statistics, and the cluster assignment for each sample.
#'
#' @examples
#'
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
#' regions <- GenomicRanges::GRanges(seqnames = c("chr19", "chr19"),
#'     ranges = IRanges::IRanges(start = c(54924104, 54874318),
#'                                 end = c(54929104, 54877536)))
#'
#' callWholeRegion(query = regions,
#'                 peaks = H3K4me3,
#'                 metadata = metadata,
#'                 normalize_columns = c("qValue", "pValue", "signalValue"),
#'                 mark = "H3K4me3",
#'                 summarize_columns = c("pValue", "qValue", "signalValue"),
#'                 heatmap = FALSE)
#'
#' @export
callWholeRegion <- function(query, peaks, metadata, mark,
                            filter = FALSE, filter_columns = NULL,
                            filter_thresholds = NULL, normalize = TRUE,
                            normalize_columns = NULL, tail = 0.005,
                            summarize_columns,
                            length = FALSE, fraction = TRUE, n = FALSE,
                            heatmap = TRUE, titles = NULL, outdir = NULL) {

    # Preprocessing
    if (filter) {

        if (is.null(filter_columns) || is.null(filter_thresholds))
            stop("Please provide names of columns to filter and specify
                thresholds to use.")

        peaks <- filterPeaks(peaks,
                            columns = filter_columns,
                            thresholds = filter_thresholds)
    }

    if (normalize) {

        if (is.null(normalize_columns))
        stop("Please provide names of columns to normalize genome-wide.")

        peaks <- normalizePeaks(peaks,
                                columns = normalize_columns,
                                tail = tail)
    }

    # Check titles
    if ((!is.null(titles)) && length(query) != length(titles))
        stop("Please provide one title per query region.")

    # Retrieve peaks: get a localPeaks object for each query region
    lpks      <- lapply(query, function(region)
                        retrievePeaks(peaks, metadata, region))

    # Construct feature matrices for each query region
    matrices  <- lapply(lpks, summarizePeaks, mark,
                        summarize_columns, length, fraction, n)

    # Convert the queries into a GRangesList in order to be able to Map over
    queries <- lapply(query, GRangesList)

    # Cluster the feature matrices
    if (!heatmap) {

        results <- Map(f = function(ft_mat, region)
            cluster(ft_mat, metadata, region, heatmap, titles, outdir),
            matrices, queries)

    } else {

        if (is.null(titles)) titles <- unlist(lapply(query, GRangesToCoord))

        results <- Map(f = function(ft_mat, region, title)
            cluster(ft_mat, metadata, region, heatmap, title, outdir),
            matrices, queries, titles)
    }

    return(data.frame(dplyr::bind_rows(results)))

}



#' callPositionAware
#'
#' One of two main functions in the \code{chromswitch} package, this function
#' detects a switch in chromatin state in one or
#' more regions given ChIP-seq peak calls for one mark, executing the entire
#' algorithm from preprocessing to evaluating the clustering results,
#' using the whole-region strategy.
#'
#' This strategy constructs a sample-by-feature matrix to use as input for
#' hierarchical clustering by first assembling the set of unique peaks observed
#' in the region across samples. Then for each unique peak, we model the
#' presence or absence of that peak in each sample, resulting in a binary
#' feature matrix.
#'
#' @param query GRanges list containing one or more genomic regions of interest
#' in which to call a switch. The output dataframe will contain one row per
#' region in \code{query}.
#' @param peaks List of GRanges objects storing peak calls for each sample,
#' where element names correspond to sample IDs
#' @param metadata A dataframe with at least two columns: "Sample" which stores
#' the sample IDs, "Condition", which stores the biological condition labels
#' of the samples
#' @param filter Optional: logical value, filter peaks based on thresholds on
#' peak statistics? Default: FALSE. The filter step is described in
#' \code{\link{filterPeaks}}.
#' @param filter_columns If \code{filter} is TRUE, a chracter vector
#' corresponding to names of columns in the peak metadata by which to filter
#' peaks. If \code{filter} is FALSE, not used.
#' @param filter_thresholds If \code{filter} is TRUE, a numeric vector
#' corresponding to lower cutoffs applied to metadata columns in order to filter
#' peaks. Provide one per column specified in \code{filter_columns}, in the same
#' order. If \code{filter} is FALSE, not used.
#' @param reduce Optional: logical value, if TRUE, reduce gaps between nearby
#' peaks in the same sample. See more at \code{\link{reducePeaks}}.
#' Default: TRUE
#' @param gap Numeric value, specifying the threshold distance for merging.
#' Peaks in the same sample which are within this many bp of each other will
#' be merged. Default: 300
#' @param p Numeric value in [0, 1] giving the fraction of reciprocal overlap
#' to require. Default: 0.4
#' @param heatmap Optional: Logical value, plot the heatmap corresponding to
#' the hierarchical clustering result? Default: TRUE
#' @param titles Optional:  if \code{heatmap} is TRUE, a character vector
#' of the same length as \code{query}, specifying the title to use when plotting
#' each heatmap (e.g. a gene name), also reused as the
#' prefix of the name of the file where the heatmap is saved. By default, the
#' title is the genomic coordinates of the region in the form "chrN:start-end"
#' @param outdir Optional: if \code{heatmap} is TRUE, the name of the directory
#' where heatmaps should be saved
#'
#' @return Data frame with one row per region in \code{query}. Contains the
#' coordinates of the region, the number of inferred clusters, the computed
#' cluster validity statistics, and the cluster assignment for each sample.
#'
#' @examples
#'
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
#' regions <- GenomicRanges::GRanges(seqnames = c("chr19", "chr19"),
#'     ranges = IRanges::IRanges(start = c(54924104, 54874318),
#'                                 end = c(54929104, 54877536)))
#'
#' callPositionAware(query = regions, peaks = H3K4me3, metadata = metadata)
#'
#' @export
callPositionAware <- function(query, peaks, metadata,
                            filter = FALSE, filter_columns = NULL,
                            filter_thresholds = NULL, reduce = TRUE,
                            gap = 300, p = 0.4,
                            heatmap = TRUE, titles = NULL, outdir = NULL) {

    # Preprocessing
    if (filter) {

        if (is.null(filter_columns) || is.null(filter_thresholds))
            stop("Please provide names of columns to filter and specify
                thresholds to use.")

        peaks <- filterPeaks(peaks,
                            columns = filter_columns,
                            thresholds = filter_thresholds)
    }

    # Retrieve peaks: get a localPeaks object for each query region
    lpks <- lapply(query, function(region)
        retrievePeaks(peaks, metadata, region))

    if (reduce) {

        lpks <- lapply(lpks, reducePeaks, gap)

    }

    # Construct feature matrices for each query region
    matrices  <- lapply(lpks, binarizePeaks, p)

    # Convert the queries into a GRangesList in order to be able to Map over
    queries <- lapply(query, GRangesList)

    # Cluster the feature matrices
    if (!heatmap) {

        results <- Map(f = function(ft_mat, region)
            cluster(ft_mat, metadata, region, heatmap, titles, outdir),
            matrices, queries)

    } else {

        if (is.null(titles)) titles <- unlist(lapply(query, GRangesToCoord))

        results <- Map(f = function(ft_mat, region, title)
            cluster(ft_mat, metadata, region, heatmap, title, outdir),
            matrices, queries, titles)
    }

    return(data.frame(dplyr::bind_rows(results)))

}
