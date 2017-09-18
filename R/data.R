#' H3K4me3 peak calls in a short region for six adult tissues
#'
#' A toy dataset containing MACS2 narrow peak calls for 3 brain tissues and
#' 3 other adult tissues from the Roadmap Epigenomics Project, restricted
#' to a short region on chromosome 19. The generation of this dataset is
#' executed by the script in the "data-raw" directory of this package, which
#' can be viewed at
#' \url{https://github.com/selinj/chromswitch/tree/master/data-raw}.
#'
#' @format A list with six entries, named according to fake IDs for the samples.
#' Each element contains a GRanges object with peak calls and associated
#' statistics which are computed by MACS2. This is the format expected by the
#' \code{peaks} argument in functions in chromswitch.
#'
#' @source \url{egg2.wustl.edu/roadmap/web_portal/}
"H3K4me3"
