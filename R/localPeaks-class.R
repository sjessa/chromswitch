# ---------------------------------------------------------------------------- #
#
# localPeaks object
# Defines an S4 class "localPeaks", accessor methods, and other associated
# methods for working with localPeaks objects
#
# ---------------------------------------------------------------------------- #


#' localPeaks
#'
#' An S4 class to aggregate and summarize peaks located in an input region
#'
#' @keywords internal
setClass("localPeaks",

        # Define slot contents
        slots = c(
            region = "GRanges",
            peaks = "list",
            samples = "character")
)


#' localPeaks object
#'
#' A localPeaks object is a container for the peaks for one or more marks for
#' a set of samples in a specific genomic region of interest, as well as the
#' genomic region itself, and the sample IDs. These components are needed to
#' convert sets of peaks into rectangular feature-by-sample matrices which we
#' can then use for downstream analysis - and in particular, as input to a
#' clustering algorithm in order to call a chromatin state switch. This function
#' is a constructor for a localPeaks object.
#'
#' @param region A GRanges object specifying one genomic region,
#' the query region
#' @param peaks List of lists of GRanges objects. Each outer list stores peaks
#' for each sample for one mark in \code{region}.
#' @param samples Character vector with sample identifiers.
#'
#' @return localPeaks object
#'
#' @keywords internal
localPeaks <- function(region, peaks, samples) {

    new("localPeaks",
        region = region,
        peaks = peaks,
        samples = samples)

}


#' lpkRegion
#'
#' Accessor for \code{region} slot of a \code{\link{localPeaks}} object.
#'
#' @param lpks localPeaks object
#'
#' @return GRanges object with query region associated with the localPeaks
#' object
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
#' lpk <- retrievePeaks(H3K4me3,
#'     metadata = metadata,
#'     region = GenomicRanges::GRanges(seqnames = "chr19",
#'     ranges = IRanges::IRanges(start = 54924104, end = 54929104)))
#'
#' lpkRegion(lpk)
#'
#' @export
lpkRegion <- function(lpks) lpks@region


#' lpkSamples
#'
#' Accessor for \code{samples} slot of a \code{\link{localPeaks}} object.
#'
#' @param lpks localPeaks object
#'
#' @return Character vector with sample IDs for the localPeaks object
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
#'     region = GenomicRanges::GRanges(seqnames = "chr19",
#'     ranges = IRanges::IRanges(start = 54924104, end = 54929104)))
#'
#' lpkSamples(lpk)
#'
#' @export
lpkSamples <- function(lpks) lpks@samples


#' lpkPeaks
#'
#' Accessor for \code{peaks} slot of a \code{\link{localPeaks}} object.
#'
#' @param lpks localPeaks object
#'
#' @return List of lists of GRanges objects. Each outer list stores peaks
#' for each sample for one mark in the region given by \code{lpkRegion(lpks)}.
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
#'     region = GenomicRanges::GRanges(seqnames = "chr19",
#'     ranges = IRanges::IRanges(start = 54924104, end = 54929104)))
#'
#' lpkPeaks(lpk)
#'
#' @export
lpkPeaks <- function(lpks) lpks@peaks


#' is.empty
#'
#' @return Logical value indicating if the object is empty (TRUE) or not (FALSE)
#'
#' @param object The object to check for emptiness
#'
#' @keywords internal
setGeneric("is.empty", function(object) {standardGeneric("is.empty")})


#' @describeIn is.empty Returns TRUE if the localPeaks object has no peaks in
#' any of the samples in \code{object@peaks}, i.e. if no peaks were found in
#' the query region.
#'
#' @keywords internal
setMethod("is.empty", signature(object = "localPeaks"),
        function(object) {sum(unlist(lapply(lpkPeaks(object), length))) == 0})
