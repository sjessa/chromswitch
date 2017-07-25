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

#' @importFrom magrittr %>%
NULL

#' @importFrom stats hclust dist cutree
NULL

#' @importFrom GenomeInfoDb seqinfo
NULL

#' @importClassesFrom GenomicRanges GRanges GRangesList
NULL
