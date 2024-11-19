#' Generate an Inverse CDF Function Using Akima Spline Interpolation
#'
#' This function creates an inverse cumulative distribution function (CDF)
#' for a given dataset using Akima spline interpolation. The resulting function
#' maps probabilities (in the range [0, 1]) to values in the dataset.
#'
#' @param X A numeric vector. The dataset for which the inverse CDF is to be created.
#' @return A function that takes a single argument, \code{p}, a numeric vector of
#'         probabilities in [0, 1], and returns the corresponding values interpolated
#'         from the dataset using Akima splines.
#' @details
#' The function works as follows:
#' \enumerate{
#'   \item Computes the empirical CDF (ECDF) of the dataset.
#'   \item Extracts the sorted ECDF values for the dataset.
#'   \item Sorts the original data values.
#'   \item Uses the \code{\link[akima]{aspline}} function to create a spline interpolation
#'         mapping probabilities to dataset values.
#' }
#' The resulting function leverages Akima splines, which are smooth and flexible for
#' interpolating data.
#'
#' @examples
#' # Example usage:
#' library(akima)
#' data <- c(1, 2, 3, 4, 5)
#' inv_cdf <- genCDFInv_akima(data)
#' inv_cdf(c(0.1, 0.5, 0.9))  # Compute interpolated values for given probabilities
#'
#' @seealso \code{\link[stats]{ecdf}}, \code{\link[akima]{aspline}}
#' @importFrom akima aspline
#' @export
genCDFInv_akima <- function(X) {
  # Compute the empirical cumulative distribution function (ECDF)
  ecdf1 <- ecdf(X)

  # Extract sorted ECDF values and data
  U <- sort(ecdf1(X))
  Finv <- sort(X)

  # Create a function using Akima spline interpolation
  ak2 <- function(x) {
    aspline(U, Finv, x)$y
  }

  return(ak2)
}
