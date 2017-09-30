# --------------------------------------------------------------------------- #
#
# Functions for dealing with genomic coordinates, simple manipulations
# of GRanges objects, etc.
#
# --------------------------------------------------------------------------- #

#' makeBrowserCoord
#'
#' Given coordinates for a genomic region, return a browser-friendly version.
#'
#' @param chr The chromosome
#' @param start The starting position of the genomic region
#' @param end The ending position of the genomic region
#'
#' @examples
#'
#' makeBrowserCoord("chr1", 1000, 2000)
#'
#' @export
#' @return String with copy-pastable, genome browser-friendly version of
#' coordinates.
makeBrowserCoord <- function(chr, start, end) {

    paste0(chr, ":", start, "-", end)

}


#' coordToGRanges
#'
#' Convert a string of genomic coordinates to a GRanges object
#'
#' @param coord String coordinate in genome browser-friendly format to convert
#' to a GRanges object
#'
#' @examples
#' string <- "chr1:1000-2000"
#' coordToGRanges(string)
#'
#' @export
#' @return GRanges object
coordToGRanges <- function(coord) {
    # String of coordinates to GRanges

    split <- strsplit(coord, ":")
    split2 <- strsplit(split[[1]][2], "-")
    GRanges(seqnames = split[[1]][1],
            ranges = IRanges(as.numeric(gsub(",", "", split2[[1]][1])),
                                    as.numeric(gsub(",", "", split2[[1]][2]))))

}


#' GRangesToCoord
#'
#' Convert a GRanges object for one region to a genome browser-friendly string
#'
#' @param gr GRanges object specifying region to convert to a string
#'
#' gr <- GRanges(seqnames = "chr1",
#'               ranges = IRanges(start = 1000, end = 2000))
#'
#' GRangesToCoord(gr)
#' @export
#' @return String
GRangesToCoord <- function(gr) {

    makeBrowserCoord(seqnames(gr), start(gr), end(gr))

}
