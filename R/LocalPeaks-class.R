# ---------------------------------------------------------------------------- #
#
# LocalPeaks object
# Defines an S4 class "LocalPeaks", accessor generics and methods,
# and other associated methods for working with LocalPeaks objects
#
# ---------------------------------------------------------------------------- #


setClass("LocalPeaks",

        # Define slot contents
        slots = c(
            region = "GRanges",
            peaks = "list",
            samples = "character")
)


# LocalPeaks object
#
# This function is a constructor for a LocalPeaks object.
#
# @param region A GRanges object specifying one genomic region,
# the query region
# @param peaks List of lists of GRanges objects. Each outer list stores peaks
# for each sample for one mark in \code{region}.
# @param samples Character vector with sample identifiers.
#
# @return LocalPeaks object
LocalPeaks <- function(region, peaks, samples) {

    new("LocalPeaks",
        region = region,
        peaks = peaks,
        samples = as.character(samples))

}


# region
#
# Generic function
#
# @param x Object
#
# @return A region associated with the object
setGeneric("region", function(x) standardGeneric("region"))


#' region
#'
#' Accessor for \code{region} slot of a \code{\linkS4class{LocalPeaks}} object.
#'
#' @param x LocalPeaks object
#'
#' @return GRanges object with query region associated with the LocalPeaks
#' object
#'
#' @examples
#'
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
#' region(lpk)
#'
#' @export
#' @aliases region-method
setMethod("region", signature(x = "LocalPeaks"),
          function(x) x@region)


# peaks
#
# Generic function
#
# @param x Object
#
# @return A set of peaks associated with the object
setGeneric("peaks", function(x) standardGeneric("peaks"))


#' peaks
#'
#' Accessor for \code{peaks} slot of a \code{\linkS4class{LocalPeaks}} object.
#'
#' @param x LocalPeaks object
#'
#' @return List of lists of GRanges objects. Each outer list stores peaks
#' for each sample for one mark in the region given by \code{region(lpks)}.
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
#' peaks(lpk)
#'
#' @export
#' @aliases peaks-method
setMethod("peaks", signature(x = "LocalPeaks"),
          function(x) x@peaks)


#' samples
#'
#' Accessor for \code{samples} slot of a \code{\linkS4class{LocalPeaks}} object.
#'
#' @param object LocalPeaks object
#'
#' @return Character vector with sample IDs for the LocalPeaks object
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
#' samples(lpk)
#'
#' @export
#' @aliases samples-method
setMethod("samples", signature(object = "LocalPeaks"),
          function(object) object@samples)


# isEmpty
#
# Returns TRUE if the LocalPeaks object has no peaks in
# any of the samples in \code{object@peaks}, i.e. if no peaks were found in
# the query region.
setMethod("isEmpty", signature(x = "LocalPeaks"),
        function(x) {sum(unlist(lapply(peaks(x), length))) == 0})
