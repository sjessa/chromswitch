# ---------------------------------------------------------------------------- #
#
# Functions related to reading and writing BED files
#
# ---------------------------------------------------------------------------- #

#' loadSampleBed
#'
#' Given a path, read a BED file into a dataframe
#'
#' @param bed_file Character string containing path to BED file with first three
#' columns "chrm", "start", "end" and optional additional columns
#' @param metadata_cols Optional: a character vector containing names of
#' all additional columns. If not specified but additional columns are found,
#' they are dropped. If some are specified, it is assumed that they describe
#' the columns which immediately follow "end".
#'
#' @return A GRanges object
#' @keywords internal
loadSampleBed <- function(bed_file, metadata_cols = NULL) {

    if (missing(metadata_cols)) columns <- c(c("chr", "start", "end"))
    else columns <- c(c("chr", "start", "end"), metadata_cols)

    bed <- readr::read_tsv(file = bed_file, col_names = columns)
    bed_gr <- makeGRangesFromDataFrame(bed, keep.extra.columns = TRUE)

    return(bed_gr)

}


#' loadBed
#'
#' Load BED files for one mark for all samples into GRanges objects.
#'
#' @param metadata A dataframe with a column "Sample" which stores
#' the sample identifiers, and at least one column, titled by the histone mark
#' or ChIP-seq target, storing paths to the BED files containing peak calls
#' @param mark String, the name of the mark for which peak calls should be
#' loaded, must correspond to a column in \code{metadata}
#' @param metadata_cols Optional: a character vector containing names of
#' all additional columns. If not specified but additional columns are found,
#' they are dropped. If some are specified, it is assumed that they describe
#' the columns which immediately follow "end".
#'
#' @return A list of GRanges objects storing peak calls for each sample.
#' List elements are named using provided sample names.
#'
#' @examples
#'
#' samples <- c("E068", "E071", "E074", "E101", "E102", "E110")
#' outfiles <- system.file("extdata", paste0(samples, ".H3K4me3.bed"),
#' package = "chromswitch")
#'
#' metadata <- data.frame(Sample = samples,
#'                         H3K4me3 = outfiles,
#'                         stringsAsFactors = FALSE)
#'
#' loadBed(metadata,
#'         mark = "H3K4me3",
#'         metadata_cols = c("name", "score","strand",
#'                             "signalValue", "pValue", "qValue", "peak"))
#'
#' @export
loadBed <- function(metadata, mark, metadata_cols = NULL) {

    if (!(mark %in% names(metadata))) stop("Paths for this mark was not found.
                                        Ensure the 'mark' argument specifies
                                        the name of a column in 'metadata'.")

    paths <- metadata[[mark]]
    beds <- lapply(paths, loadSampleBed, metadata_cols)
    names(beds) <- metadata$Sample

    return(beds)

}

