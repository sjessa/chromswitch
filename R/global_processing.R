# ---------------------------------------------------------------------------- #
#
# Global processing functions:
# Functions which take in metadata or global (i.e. genome-wide) peak calls and
# apply manipulations globally, returning global data
#
# ---------------------------------------------------------------------------- #


#' winsorNorm
#'
#' Normalize a numeric vector by rescaling and Winsorizing, i.e. rescale the
#' middle of the data to the range [0, 1] and bound the upper tail to 1 and the
#' lower tail to 0, effectively replacing a fixed amount of extreme values in
#' each tail. Similar to trimming the tails except instead of discarding
#' the tails entirely they're bounded.
#'
#' @param x A numeric vector, the data to be normalized
#' @param trim Numeric, a fraction in [0, 1] specifying how much of the data
#' to bound to 0 (for the lower tail) or 1 (for the upper tail)
#'
#' @return Numeric vector
#' @keywords internal
winsorNorm <- function(x, trim) {

    trimmed    <- DescTools::Trim(x, trim = trim)
    scaled     <- (x - min(trimmed))/diff(range(trimmed))
    winsorized <- DescTools::Winsorize(scaled, minval = 0, maxval = 1)

    return(winsorized)

}


#' normalizePeaks
#'
#' Given a set of peak calls for different marks and samples as returned by
#' \code{\link{loadBed}},, normalize all peaks genome-wide for each
#' sample and mark by rescaling and Winsorizing,
#' i.e. rescale the middle of the data to the range
#' [0, 1] and bound the upper tail to 1 and the lower tail to 0, effectively
#' replacing a fixed amount of extreme values in each tail. Similar to trimming
#' the tails except instead of discarding the tails entirely they're bounded.
#'
#' @param peaks List of lists of GRanges objects, as returned by
#' \code{\link{loadBed}}.
#' @param columns Character vector specifying the names of columns to normalize
#' @param tail Optional: numeric, a fraction in [0, 1] specifying how much of
#' the data to bound to 0 (for the lower tail) or 1 (for the upper tail).
#' Default: 0.005.
#'
#' @examples
#' normalizePeaks(H3K4me3, columns = c("signalValue", "pValue", "qValue"))
#'
#' @return A list of GRanges objects storing peak calls for each sample, with
#' columns specified in \code{columns} normalized.
#'
#' @export
#' @seealso winsorNorm
normalizePeaks <- function(peaks, columns, tail = 0.005) {

    if (!all(columns %in% names(mcols(peaks[[1]])))) stop("One or more columns
                                    was not found in the metadata columns for
                                    the provided peaks.")

    if (tail < 0 || tail > 1) stop("Please provide a value for 'tail'
                                    between 0 and 1.")

    # A function to normalize all specified columns in a dataframe
    normalizeStats <- function(df) dplyr::mutate_at(df,
                                            .vars = columns,
                                            trim = tail,
                                            .funs = winsorNorm)

    peaks %>%
        lapply(as.data.frame) %>%
        lapply(normalizeStats) %>%
        lapply(makeGRangesFromDataFrame, keep.extra.columns = TRUE)

}


#' filterSamplePeaks
#'
#' @param sample_peaks GRanges object (peak calls for one sample)
#' @param columns Character vector of column names containing stats by which to
#' filter peaks
#' @param thresholds Vector of numeric values giving the lower thresholds to use
#' for each of the columns specified, in the same order as \code{columns}
#'
#' @return A GRanges object containing peak calls for this sample, filtered
#' according to the thresholds and columns specified.
#'
#' @seealso filterPeaks
#' @keywords internal
filterSamplePeaks <- function(sample_peaks, columns, thresholds) {

    pk <- as.data.frame(sample_peaks)

    for (i in 1:length(columns)) {

        # NSE workaround using lazyeval::interp
        criteria <- lazyeval::interp(~ stat >= thresholds[i],
                                    stat = as.name(columns[i]))
        pk <- pk %>% dplyr::filter_(criteria)

    }

    return(makeGRangesFromDataFrame(pk, keep.extra.columns = TRUE))
}


#' filterPeaks
#'
#' Given a set of peak calls for different marks and samples as returned by
#' \code{\link{loadBed}}, filter peaks according to values in numeric
#'
#' @param peaks List of GRanges objects, as returned by \code{\link{loadBed}}
#' @param columns Character vector of column names containing stats by which to
#' filter peaks
#' @param thresholds Vector of numeric values giving the lower thresholds to use
#' for each of the columns specified, in the same order as \code{columns}
#'
#' @return A list of GRanges objects storing peak calls for each sample, with
#' peaks filtered according to the columns and thresholds specified.
#'
#' @examples
#' filterPeaks(peaks = H3K4me3,
#'     columns = c("signalValue", "pValue"),
#'     thresholds = c(4, 10))
#'
#' @export
filterPeaks <- function(peaks, columns, thresholds) {

    if (length(columns) != length(thresholds))
    stop("Please provide one threshold per column in order to filter peaks.")

    if (!all(columns %in% names(mcols(peaks[[1]])))) stop("One or more columns
                                    was not found in the metadata columns for
                                    the provided peaks.")

    cols_are_numeric <- mcols(peaks[[1]]) %>%
        as.data.frame %>%
        dplyr::select_(.dots = columns) %>%
        lapply(is.numeric) %>%
        unlist()

    if (!all(cols_are_numeric)) stop("One or more columns specified in 'columns'
                            is not numeric. Only numeric columns can
                            be filtered.")

    peaks %>% lapply(filterSamplePeaks, columns, thresholds)

}


