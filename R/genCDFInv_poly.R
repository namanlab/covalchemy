#' Generate an Inverse CDF Function Using Polynomial Regression
#'
#' This function creates an inverse cumulative distribution function (CDF) for a
#' given dataset using polynomial regression. The resulting function maps probabilities
#' (in the range [0, 1]) to values in the dataset.
#'
#' @param data A numeric vector. The dataset for which the inverse CDF is to be created.
#' @param degree An integer. The degree of the polynomial to fit in the regression.
#' @return A function that takes a single argument, \code{y}, a numeric vector of
#'         probabilities in [0, 1], and returns the corresponding values predicted
#'         by the polynomial regression.
#' @details
#' The function works as follows:
#' \enumerate{
#'   \item Sorts the dataset and computes the empirical CDF (ECDF) of the data.
#'   \item Fits a polynomial regression to model the relationship between
#'         ECDF values and the sorted dataset.
#'   \item Uses the fitted polynomial model to predict the inverse CDF for
#'         given probabilities.
#' }
#' The degree of the polynomial can be specified to control the flexibility
#' of the regression model.
#'
#' @examples
#' # Example usage:
#' data <- c(1, 2, 3, 4, 5)
#' inv_cdf <- genCDFInv_poly(data, degree = 2)
#' inv_cdf(c(0.1, 0.5, 0.9))  # Compute predicted values for given probabilities
#'
#' @seealso \code{\link[stats]{lm}}, \code{\link[stats]{poly}}, \code{\link[stats]{ecdf}}
#' @export
genCDFInv_poly <- function(data, degree) {
  # Sort the data and compute the empirical cumulative distribution function (ECDF)
  x <- sort(data)
  y <- ecdf(data)(x)

  # Fit polynomial regression model of specified degree
  poly_fit <- lm(x ~ poly(y, degree = degree, raw = TRUE))

  # Create a function to predict values for given probabilities using the polynomial model
  cdf_poly <- function(y) predict(poly_fit, newdata = data.frame(y = y))

  return(cdf_poly)
}
