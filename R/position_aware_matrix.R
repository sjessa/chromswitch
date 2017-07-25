# ---------------------------------------------------------------------------- #
#
# Feature matrix construction functions for the position-aware strategy:
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
#' @return Logical value, TRUE if \code{a} and \code{b} are the same
#' by having a p-reciprocal overlap, FALSE otherwise
#'
#' @examples
#' a <- GenomicRanges::GRanges(seqnames = "chr11",
#'             ranges = IRanges::IRanges(start = 112829468, end = 112834468))
#'
#' b <- GenomicRanges::GRanges(seqnames = "chr11",
#'             ranges = IRanges::IRanges(start = 112829468, end = 113834468))
#'
#' c <- GenomicRanges::GRanges(seqnames = "chr11",
#'             ranges = IRanges::IRanges(start = 112829968, end = 112834968))
#'
#' # Expect FALSE
#' pReciprocalOverlap(a, b, 0.9)
#'
#' # Expect TRUE
#' pReciprocalOverlap(a, c, 0.9)
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


#' getUniquePeaks
#'
#' Given a set of peaks, collapse them such that we retain only unique peaks
#' which do not have a reciprocal overlap of \code{p}%. We do so by scanning
#' through the union of peaks across samples in the query region in the order
#' they're provided, and growing a list of unique peaks. For every peak visited,
#' add it to the list of unique peaks if it doesn't have a p-reciprocal overlap
#' with any peak already in the list;E otherwise, discard it.
#'
#' @param loc_union GRanges object storing union of peaks in the query
#' region across samples
#' @param p Numeric value in [0, 1] giving the fraction of reciprocal overlap
#' to require.
#'
#' @return GRanges objects with a list of unique peaks as determined by
#' the p-reciprocal overlap rule
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



