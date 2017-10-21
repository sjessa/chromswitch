# ---------------------------------------------------------------------------- #
#
# Feature matrix construction functions for the summary strategy:
# Functions which take in a LocalPeaks object or a sample-by-feature matrix,
# and return a feature matrix, and their helpers
#
# ---------------------------------------------------------------------------- #


# peakOverlap
#
# Compute the fraction of the query region of interest which is overlapped
# by peaks. This value will be at most 1, since we only consider the
# overlapping portion of each peak, and peaks do not overlap with each other.
# This value will be 0 in the case that there are no peaks in the query region.
#
# @param region GRanges object defining query region
# @param peaks GRanges object storing peak calls for a sample in the
# query region
#
# @return double, a number in [0,1].
peakOverlap <- function(region, peaks) {

    # If no peaks found to be overlapping, return 0
    if (length(peaks) == 0) {

        return(0)

    } else {
        # If a peak is only partially overlapping the query
        # region, trim its start/endpoints to match the region exactly
        contained_regions <- peaks %>%
            as.data.frame %>%
            dplyr::mutate(start = ifelse(start < start(region),
                                        start(region),
                                        start)) %>%
            dplyr::mutate(end = ifelse(end > end(region),
                                        end(region),
                                        end)) %>%
            makeGRangesFromDataFrame(keep.extra.columns = TRUE)

        fraction <- sum(width(contained_regions))/width(region)
        return(fraction)
    }
}


# summarizeSamplePeaks
#
# For one sample, given peaks in the query region, compute several summary
# statistics of peaks in that region, including the mean, median, and max of
# any specified
#
# @param peaks GRanges object storing peak calls for a sample in the
# query region
# @param query GRanges object storing query region
# @param mark String specifying the name of the mark for which peaks are given
# @param cols Character vector of column names on which to compute summary
# statistics
# @param fraction Loogical: compute the fraction of the region overlapped by
# peaks?
# @param n Logical: compute the number of peaks in the region?
#
# @return data.frame
summarizeSamplePeaks <- function(peaks, query, mark, cols =  NULL,
                                fraction = TRUE, n = FALSE) {

    cols_are_numeric <- mcols(peaks) %>%
        as.data.frame %>%
        dplyr::select_(.dots = cols) %>%
        lapply(is.numeric) %>%
        unlist()

    if (!all(cols_are_numeric)) stop("One or more columns specified in 'cols' is
                            not numeric. Summary statistics can only be computed
                            on numeric data.")

    # Compute statistics summarizing the peaks located in the query region
    # for specified columns & length
    if (length(peaks) > 0) {

        stats <- mcols(peaks) %>%
                as.data.frame %>%
                dplyr::summarise_at(.vars = c("length", cols),
                                    # Summary stats to calculate
                                    dplyr::funs(mean, median, max)) %>%
                # Compute fraction of the region which is overlapped by peaks
                dplyr::mutate(fraction_region_in_peaks = peakOverlap(query,
                                                                     peaks)) %>%
                # Compute number of peaks in region
                dplyr::mutate(n_peaks = length(peaks))

    } else if (length(peaks) == 0) {

        stats <- rep(0, length(c("length", cols)))
        names(stats) <- c("length", cols)

        # No peaks in the query region, set all values to 0
        stats <- stats %>%
            t() %>%
            magrittr::set_rownames(NULL) %>%
            as.data.frame %>%
            dplyr::mutate_all(dplyr::funs(mean, median, max)) %>%
            dplyr::select_(.dots = paste("-", c("length", cols))) %>%
            dplyr::mutate(fraction_region_in_peaks = peakOverlap(query,
                                                                peaks)) %>%
            dplyr::mutate(n_peaks = length(peaks))
    }

    # Optional additional summary statistics
    if (!n)        stats <- stats %>% dplyr::select_("-n_peaks")
    if (!fraction) stats <- stats %>%
            dplyr::select_("-fraction_region_in_peaks")

    stats <- stats %>% dplyr::select(-dplyr::matches("length"))
    if(is.null(cols)) stats <- stats %>%
        dplyr::select(-dplyr::one_of(c("mean", "median", "max")))


    names(stats) <- paste0(mark, "_", names(stats))

    return(stats)

}


#' summarizePeaks
#'
#' Given peaks for a set of samples in a query region, construct a sample-by-
#' feature matrix where each row is a vector of summary statistics computed from
#' peaks in the region.
#'
#' @param localpeaks LocalPeaks object
#' @param mark String specifying the name of the mark for which the LocalPeaks
#' object is given
#' @param cols Character vector of column names on which to compute summary
#' statistics
#' @param fraction Loogical: compute the fraction of the region overlapped by
#' peaks?
#' @param n Logical: compute the number of peaks in the region?
#'
#' @return A matrix where rows are samples and columns are features
#'
#' @examples
#' samples <- c("E068", "E071", "E074", "E101", "E102", "E110")
#' bedfiles <- system.file("extdata", paste0(samples, ".H3K4me3.bed"),
#' package = "chromswitch")
#'
#' metadata <- data.frame(Sample = samples,
#'     H3K4me3 = bedfiles,
#'     stringsAsFactors = FALSE)
#'
#' lpk <- retrievePeaks(H3K4me3,
#'     metadata = metadata,
#'     region = GRanges(seqnames = "chr19",
#'     ranges = IRanges(start = 54924104, end = 54929104)))
#'
#' summarizePeaks(lpk, mark = "H3K4me3", cols = c("qValue", "signalValue"))
#'
#' @export
summarizePeaks <- function(localpeaks, mark, cols,
                            fraction = TRUE, n = FALSE) {

    if(is.null(cols) && fraction == FALSE && n == FALSE) {
        stop("No columns specified, and none of length, fraction, or
             n set to TRUE, therefore cannot construct a feature matrix.")
    }

    ft_matrix <- peaks(localpeaks) %>%
        lapply(summarizeSamplePeaks,
                region(localpeaks), mark, cols, fraction, n) %>%
        dplyr::bind_rows()

    rownames(ft_matrix) <- samples(localpeaks)
    return(data.matrix(ft_matrix))

}

