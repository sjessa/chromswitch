# ---------------------------------------------------------------------------- #
#
# Helper functions for reading in files
#
# ---------------------------------------------------------------------------- #

#' readNarrowPeak
#'
#' A helper function for reading in narrow peak calls for a set of samples.
#' Peak calls are assumed to be in ENCODE narrowPeak format
#' (\url{https://genome.ucsc.edu/FAQ/FAQformat.html#format12}) as returned by
#' MACS2 (\url{http://liulab.dfci.harvard.edu/MACS/}). This is BED6+4 format.
#'
#' @param paths Character vector storing paths for BED files containing peak
#' calls for each sample, in the same order as in the Sample column of
#' \code{metadata}.
#' @param metadata A dataframe with at least two columns: "Sample" which stores
#' the sample identifiers, and "Condition" which stores
#' the biological condition labels of the samples.
#'
#' @return Named list of GRanges objects containing peak calls for each sample.
#' @export
#'
#' @examples
#'
#' samples <- c("E068", "E071", "E074", "E101", "E102", "E110")
#' bedfiles <- system.file("extdata", paste0(samples, ".H3K4me3.bed"),
#' package = "chromswitch")
#' Conditions <- c(rep("Brain", 3), rep("Other", 3))
#'
#' metadata <- data.frame(Sample = samples,
#'     H3K4me3 = bedfiles,
#'     Condition = Conditions,
#'     stringsAsFactors = FALSE)
#'
#' readNarrowPeak(bedfiles, metadata)
readNarrowPeak <- function(paths, metadata) {

    extra_cols <- c("name" = "character",
                    "score" = "integer",
                    "strand" = "character",
                    "signalValue" = "numeric",
                    "pValue" = "numeric",
                    "qValue" = "numeric",
                    "peak" = "numeric")

    pk <- lapply(paths, rtracklayer::import,
                 format = "bed", extraCols = extra_cols)
    names(pk) <- metadata$Sample

    return(pk)

}
