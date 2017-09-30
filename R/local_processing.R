# ---------------------------------------------------------------------------- #
#
# Local processing functions:
# Functions which take in global sets of peak calls or localPeaks objects,
# performing some manipulation on the local peaks, returning localPeaks objects
#
# ---------------------------------------------------------------------------- #


# retrieveSamplePeaks
#
# Given a set of peaks for one sample and a region of interest,
# retrieve peaks overlapping query region
#
# @param peaks GRanges object containing peak calls for the sample
# @param region GRanges object specifying genomic region in which to search
# for peaks
#
# @return GRanges object. Contains all peaks overlapping query region.
retrieveSamplePeaks <- function(peaks, region) {

    # Retrieve peaks and associated metadata columns in query region
    olaps       <- as.data.frame(GenomicRanges::findOverlaps(region, peaks))
    idx         <- olaps$subjectHits
    local_peaks <- peaks[idx, ] # Not a localPeaks object, just the peaks

    # Add peak length as a metadata column using the GR width accessor
    mcols(local_peaks) <- c(mcols(local_peaks),
                            data.frame(length = width(local_peaks)))

    return(local_peaks)

}


#' retrievePeaks
#'
#' Given a peak calls for a set of samples, for each sample, get the peaks which
#' overlap a specified genomic region of interest. Typically, this corresponds
#' to the region for which we will construct a feature matrix representing
#' peaks in the region in order to call a chromatin state switch.
#'
#' @param peaks List of GRanges objects storing peak calls for each sample
#' @param metadata Dataframe with a column "Sample" which stores
#' the sample identifiers, and at least one column, titled by the histone mark
#' or ChIP-seq target, storing paths to the BED files containing peak calls
#' @param region GRanges object specifying one genomic region,
#' the query region
#'
#' @return localPeaks object as described in \code{\linkS4class{localPeaks}}
#'
#' @examples
#'
#' samples <- c("E068", "E071", "E074", "E101", "E102", "E110")
#' outfiles <- system.file("extdata", paste0(samples, ".H3K4me3.bed"),
#' package = "chromswitch")
#'
#' metadata <- data.frame(Sample = samples,
#'     H3K4me3 = outfiles,
#'     stringsAsFactors = FALSE)
#'
#' retrievePeaks(H3K4me3,
#'     metadata = metadata,
#'     region = GRanges(seqnames = "chr19",
#'     ranges = IRanges(start = 54924104, end = 54929104)))
#'
#' @export
retrievePeaks <- function(peaks, metadata, region) {

    local_peaks <- lapply(peaks, retrieveSamplePeaks, region)
    localPeaks(region  = region,
                peaks   = local_peaks,
                samples = metadata$Sample)

}


#' reducePeaks
#'
#' Given a localPeaks object, merge peaks which are in the same sample and are
#' separated by no more than \code{gap} base pairs. When two non-overlapping
#' peaks are merged, a new peak is created which starts at the starting position
#' of the first peak and ends at the ending position of the second peak,
#' spanning the range of both peaks and the gap between them.
#'
#' @param localpeaks localPeaks object
#' @param gap Numeric value, specifying the threshold distance for merging.
#' Peaks in the same sample which are within this many bp of each other will
#' be merged.
#'
#' @return The localPeaks object that was provided as input, with nearby peaks
#' merged
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
#' reducePeaks(lpk, gap = 300)
#'
#' @export
reducePeaks <- function(localpeaks, gap) {

    if (gap <= 0) stop("The gap argument must be a positive integer.")

    localpeaks@peaks <- peaks(localpeaks) %>% lapply(GenomicRanges::reduce,
                                        min.gapwidth = gap + 1)
    return(localpeaks)

}

