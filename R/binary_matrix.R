# ---------------------------------------------------------------------------- #
#
# Feature matrix construction functions for the binary strategy:
# Functions which take in a localPeaks object or a sample-by-feature matrix,
# and return a feature matrix, and their helpers
#
# ---------------------------------------------------------------------------- #


#' pReciprocalOverlap
#'
#' If a and b denote two genomic regions, check whether they overlap
#' reciprocally by \code{p}*100% of their lengths
#'
#' @param a GRanges object storing first region
#' @param b GRanges object storing second region
#' @param p Numeric value in [0, 1] giving the fraction of reciprocal overlap
#' to require.
#'
#' @examples
#' a <- GRanges(seqnames = "chr11",
#'              ranges = IRanges(start = 112829468, end = 112834468))
#' b <- GRanges(seqnames = "chr11",
#'              ranges = IRanges(start = 112829468, end = 113834468))
#'
#' pReciprocalOverlap(a, b, 0.9)
#'
#' @return Logical value, TRUE if \code{a} and \code{b} are the same
#' by having a p-reciprocal overlap, FALSE otherwise
#'
#' @export
pReciprocalOverlap <- function(a, b, p) {

    isOlap <- function(b_i, a, p) {

        # Require the larger amount of overlap, effectively a p-reciprocal olap
        minolap <- max(p * width(a), p * width(b_i))
        IRanges::overlapsAny(a, b_i, minoverlap = minolap)

    }

    any(unlist(lapply(b, isOlap, a, p)))

}


# getUniquePeaks
#
# Given a set of peaks, collapse them such that we retain only unique peaks
# which do not have a reciprocal overlap of \code{p}%. We do so by scanning
# through the union of peaks across samples in the query region in the order
# they're provided, and growing a list of unique peaks. For every peak visited,
# add it to the list of unique peaks if it doesn't have a p-reciprocal overlap
# with any peak already in the list;E otherwise, discard it.
#
# @param loc_union GRanges object storing union of peaks in the query
# region across samples
# @param p Numeric value in [0, 1] giving the fraction of reciprocal overlap
# to require.
#
# @return GRanges objects with a list of unique peaks as determined by
# the p-reciprocal overlap rule
getUniquePeaks <- function(loc_union, p) {

    unique_peaks <- GRanges()

    for (i in 1:length(loc_union)) {

        # Check if the current peak has a p reciprocal overlap with
        # any peaks in the list of unique peaks
        pOlap <- lapply(unique_peaks, pReciprocalOverlap, loc_union[i], p) %>%
            unlist() %>%
            any()

        # If not, add it to the list
        if (pOlap) next
        else unique_peaks <- append(unique_peaks, loc_union[i])

    }

    return(unique_peaks)

}


# getSamplePeakProfile
#
# For a set of peaks in one sample and a set of windows, which could correspond
# to bins of a region, or unique peaks, model presence or absence of peaks in
# each window by calling peaks present in a window if peaks in the sample have
# a p-reciprocal overlap with the window. Alternatively, the windows can be
# viewed as unique peaks, and this function calls the presence or absence of
# each peak in the sample.
#
# @param peaks GRanges object containing peaks for one sample
# @param windows GRanges object containing windows in which to model presence/
# absence of peaks, these become the features or columns of the output
# or columns
# @param p Numeric value in [0, 1] giving the fraction of reciprocal overlap
# to require.
#
# @return A one-row dataframe with a logical value TRUE/FALSE in each column
# (window) indicating whether any peaks overlap the window
getSamplePeakProfile <- function(peaks, windows, p) {

    overlaps <- lapply(windows, pReciprocalOverlap, peaks, p) %>%
        unlist() %>%
        t() %>%
        data.frame()

    return(overlaps)

}


#' binarizePeaks
#'
#' Given peaks for a set of samples in a query region, construct a sample-by-
#' feature matrix where each row is a binary vector which models the presence
#' or absence of unqiue peaks in the region.
#'
#' @param localpeaks localPeaks object storing peaks for all samples in the
#' query region
#' @param p Numeric value in [0, 1] giving the fraction of reciprocal overlap
#' to require.
#'
#' @examples
#' samples <- c("E068", "E071", "E074", "E101", "E102", "E110")
#' outfiles <- system.file("extdata", paste0(samples, ".H3K4me3.bed"),
#' package = "chromswitch")
#'
#' metadata <- data.frame(Sample = samples,
#'     H3K4me3 = outfiles,
#'     stringsAsFactors = FALSE)
#'
#' lpk <- retrievePeaks(H3K4me3,
#'     metadata = metadata,
#'     region = GRanges(seqnames = "chr19",
#'     ranges = IRanges(start = 54924104, end = 54929104)))
#'
#' # Get feature matrix
#' ft_matrix <- binarizePeaks(lpk, 0.5)
#'
#' # See features
#' attr(ft_matrix, "features")
#'
#' @return A data frame where rows are samples and columns are features. The
#' genomic ranges which give the features are returned as the \code{features}
#' attribute of the data frame.
#'
#' @export
binarizePeaks <- function(localpeaks, p) {

    # Empty peaks scenario: make a dummy feature, set to TRUE for all samples
    if (isEmpty(localpeaks)) {

        warning("No peaks found in region")

        ft_matrix <- data.frame(no_peak = rep(1, length(peaks(localpeaks))))
        rownames(ft_matrix) <- samples(localpeaks)

        return(ft_matrix)
    }

    # Get the union of all peaks
    loc_union <- Reduce("c", peaks(localpeaks))

    # From the union, extract unique peaks to serve as features
    uniq_pks <- getUniquePeaks(loc_union, p)

    # Model presence/absence of each feature peak in each sample
    ft_matrix <- lapply(peaks(localpeaks), getSamplePeakProfile, uniq_pks, p) %>%
        dplyr::bind_rows()

    # Convert from logical to numeric
    ft_matrix <- ft_matrix * 1

    uniq_pks_coords <- uniq_pks %>% lapply(GRangesToCoord) %>% unlist()
    colnames(ft_matrix) <- uniq_pks_coords
    rownames(ft_matrix) <- samples(localpeaks)

    # This is not necessary if the features are unique peaks rather than
    # arbitrary bins
    # # Remove columns (features) which are all FALSE i.e. no pks in region
    # ft_matrix <- as.matrix(ft_matrix) * 1
    # ft_matrix <- ft_matrix[, which(!apply(ft_matrix, 2,
    #                                       FUN = function(x){all(x == 0)}))]

    # Return the feature regions as an attribute
    attr(ft_matrix, "features") <- uniq_pks

    return(ft_matrix)

}
