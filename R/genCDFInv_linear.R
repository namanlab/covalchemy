#' Generate an Inverse CDF Function Using Linear Interpolation
#'
#' This function creates an inverse cumulative distribution function (CDF) for
#' a given dataset using linear interpolation. The resulting function maps
#' probabilities (in the range [0, 1]) to values in the dataset.
#'
#' @param X A numeric vector. The dataset for which the inverse CDF is to be created.
#' @return A function that takes a single argument, \code{p}, a numeric vector of
#'         probabilities in [0, 1], and returns the corresponding values interpolated
#'         from the dataset.
#' @details
#' The function works as follows:
#' \enumerate{
#'   \item Computes the empirical CDF (ECDF) of the dataset.
#'   \item Extracts the sorted ECDF values for the dataset.
#'   \item Sorts the original data values.
#'   \item Uses \code{\link[stats]{approxfun}} to create a linear interpolation function
#'         mapping probabilities to dataset values.
#' }
#' The resulting function can handle probabilities outside [0, 1] using the
#' \code{rule = 2} parameter in \code{\link[stats]{approxfun}}, which extrapolates
#' based on the nearest data points.
#'
#' @examples
#' # Example usage:
#' data <- c(1, 2, 3, 4, 5)
#' inv_cdf <- genCDFInv_linear(data)
#' inv_cdf(c(0.1, 0.5, 0.9))  # Compute the interpolated values for given probabilities
#'
#' @seealso \code{\link[stats]{ecdf}}, \code{\link[stats]{approxfun}}
#' @export
genCDFInv_linear <- function(X) {
  # Compute the empirical cumulative distribution function (ECDF)
  ecdf1 <- ecdf(X)

  # Extract sorted ECDF values and data
  U <- sort(ecdf1(X))
  Finv <- sort(X)

  # Create a linear interpolation function
  ak2 <- approxfun(U, Finv, method = "linear", rule = 2)

  return(ak2)
}
