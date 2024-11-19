#' Generate an Inverse CDF Function Using Quantiles
#'
#' This function creates an inverse cumulative distribution function (CDF)
#' for a given dataset using quantiles. The resulting function maps probabilities
#' (in the range [0, 1]) to quantiles of the data.
#'
#' @param data A numeric vector. The dataset for which the inverse CDF is to be created.
#' @param type Integer. Specifies the algorithm used to compute quantiles.
#'             See \code{\link[stats]{quantile}} for details. Default is \code{1}.
#' @return A function that takes a single argument, \code{p}, a numeric vector of
#'         probabilities in [0, 1], and returns the corresponding quantiles from
#'         the data.
#' @details
#' The function works by wrapping the \code{\link[stats]{quantile}} function. The
#' \code{type} parameter controls the quantile computation method. For example:
#' \itemize{
#'   \item Type 1 corresponds to inverse of the empirical distribution function (default).
#'   \item Other types correspond to different quantile algorithms as documented in
#'         \code{\link[stats]{quantile}}.
#' }
#'
#' @examples
#' # Example usage:
#' data <- c(1, 2, 3, 4, 5)
#' inv_cdf <- genCDFInv_quantile(data, type = 1)
#' inv_cdf(c(0.25, 0.5, 0.75))  # Compute the 25th, 50th, and 75th percentiles
#'
#' @seealso \code{\link[stats]{quantile}}
#' @export
genCDFInv_quantile <- function(data, type = 1) {
  # Return a function that computes quantiles based on the input probabilities
  return(function(p) {
    # Validate that probabilities are in [0, 1]
    if (any(p < 0 | p > 1)) {
      stop("All probabilities 'p' must be in the range [0, 1].")
    }
    quantile(data, p, type = type)
  })
}
