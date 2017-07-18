# --------------------------------------------------------------------------- #
#
# Functions for dealing with genomic coordinates, simple manipulations
# of GRanges objects, etc.
#
# --------------------------------------------------------------------------- #

#' makeFriendlyCoord
#'
#' Given coordinates for a genomic region, return a browser-friendly version.
#'
#' @param chr The chromosome
#' @param start The starting position of the genomic region
#' @param end The ending position of the genomic region
#'
#' @return String with copy-pastable, genome browser-friendly version of
#' coordinates.
makeFriendlyCoord <- function(chr, start, end) {

    paste0(chr, ":", start, "-", end)

}
