# ---------------------------------------------------------------------------- #
#
# Local processing functions:
# Functions which take in global sets of peak calls or localPeaks objects,
# performing some manipulation on the local peaks, returning localPeaks objects
#
# ---------------------------------------------------------------------------- #


#' retrieveSamplePeaks
#'
#' Given a set of peaks for one sample and a region of interest,
#' retrieve peaks overlapping query region
#'
#' @param peaks GRanges object containing peak calls for the sample
#' @param region GRanges object specifying genomic region in which to search
#' for peaks
#'
#' @return GRanges object. Contains all peaks overlapping query region.
retrieveSamplePeaks <- function(peaks, region) {

    # Retrieve peaks and associated metadata columns in query region
    olaps       <- as.data.frame(GenomicRanges::findOverlaps(region, peaks))
    idx         <- olaps$subjectHits
    local_peaks <- peaks[idx, ] # Not a localPeaks object, just the peaks

    # Add peak length as a metadata column using the GR width accessor for
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
#' @return localPeaks object as described in \code{\link{localPeaks}}
#'
#' @examples
#'
#' samples <- c("brain1", "brain2", "brain3", "other1", "other2", "other3")
#' outfiles <- system.file("extdata", paste0(samples, ".H3K4me3.bed"),
#' package = "chromswitch")
#'
#' metadata <- data.frame(Sample = samples,
#'     H3K4me3 = outfiles,
#'     stringsAsFactors = FALSE)
#'
#' retrievePeaks(H3K4me3,
#'     metadata = metadata,
#'     region = GenomicRanges::GRanges(seqnames = "chr19",
#'     ranges = IRanges::IRanges(start = 54924104, end = 54929104)))
#'
#' @export
retrievePeaks <- function(peaks, metadata, region) {

    local_peaks <- lapply(peaks, retrieveSamplePeaks, region)
    localPeaks(region  = region,
                peaks   = local_peaks,
                samples = metadata$Sample)

}