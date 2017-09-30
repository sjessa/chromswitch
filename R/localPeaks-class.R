# ---------------------------------------------------------------------------- #
#
# localPeaks object
# Defines an S4 class "localPeaks", accessor generics and methods,
# and other associated methods for working with localPeaks objects
#
# ---------------------------------------------------------------------------- #


setClass("localPeaks",

        # Define slot contents
        slots = c(
            region = "GRanges",
            peaks = "list",
            samples = "character")
)


# localPeaks object
#
# This function is a constructor for a localPeaks object.
#
# @param region A GRanges object specifying one genomic region,
# the query region
# @param peaks List of lists of GRanges objects. Each outer list stores peaks
# for each sample for one mark in \code{region}.
# @param samples Character vector with sample identifiers.
#
# @return localPeaks object
localPeaks <- function(region, peaks, samples) {

    new("localPeaks",
        region = region,
        peaks = peaks,
        samples = samples)

}


#' region
#'
#' Generic function
#'
#' @param x Object
#'
#' @return A region associated with the object
setGeneric("region", function(x) {standardGeneric("region")})


#' region
#'
#' Accessor for \code{region} slot of a \code{\linkS4class{localPeaks}} object.
#'
#' @param x localPeaks object
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
#'     region = GRanges(seqnames = "chr19",
#'     ranges = IRanges(start = 54924104, end = 54929104)))
#'
#' region(lpk)
#'
#' @export
#' @aliases region-method
setMethod("region", signature(x = "localPeaks"),
          function(x) x@region)


#' peaks
#'
#' Generic function
#'
#' @param x Object
#'
#' @return A set of peaks associated with the object
setGeneric("peaks", function(x) {standardGeneric("peaks")})


#' peaks
#'
#' Accessor for \code{peaks} slot of a \code{\linkS4class{localPeaks}} object.
#'
#' @param x localPeaks object
#'
#' @return List of lists of GRanges objects. Each outer list stores peaks
#' for each sample for one mark in the region given by \code{region(lpks)}.
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
#' peaks(lpk)
#'
#' @export
#' @aliases peaks-method
setMethod("peaks", signature(x = "localPeaks"),
          function(x) x@peaks)


#' samples
#'
#' Accessor for \code{samples} slot of a \code{\linkS4class{localPeaks}} object.
#'
#' @param object localPeaks object
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
#'     region = GRanges(seqnames = "chr19",
#'     ranges = IRanges(start = 54924104, end = 54929104)))
#'
#' samples(lpk)
#'
#' @export
#' @aliases samples-method
setMethod("samples", signature(object = "localPeaks"),
          function(object) object@samples)


# isEmpty
#
# Returns TRUE if the localPeaks object has no peaks in
# any of the samples in \code{object@peaks}, i.e. if no peaks were found in
# the query region.
setMethod("isEmpty", signature(x = "localPeaks"),
        function(x) {sum(unlist(lapply(peaks(x), length))) == 0})
