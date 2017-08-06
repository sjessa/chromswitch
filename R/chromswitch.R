# ---------------------------------------------------------------------------- #
#
# Package documentation and package-level imports
#
# ---------------------------------------------------------------------------- #


#' chromswitch: An R package for detecting chromatin state switches
#'
#' chromswitch implements a flexible method to detect chromatin state
#' switches between samples in two biological conditions in a specific genomic
#' region of interest given peaks called from ChIP-seq data.
#'
#' @docType package
#' @name chromswitch
NULL

#' @import methods GenomicRanges IRanges
NULL

#' @importFrom BiocParallel bplapply bpmapply bpparam
NULL

#' @importFrom magrittr %>%
NULL

#' @importFrom stats hclust dist cutree
NULL

#' @importClassesFrom GenomicRanges GRanges GRangesList
NULL

## Handle R CMD check NOTE re: the .'s that appear in pipelines
## and other undefined variables in tidyr/dplyr functions
## as per https://stackoverflow.com/a/12429344
if(getRversion() >= "2.15.1") utils::globalVariables(c(".", "Var1", "Var2",
                                                        "Freq", "Cluster",
                                                        "Condition", "Sample",
                                                        "C1", "C2"))
