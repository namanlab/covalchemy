#' Generate Gaussian Copula Samples for Two Variables
#'
#' This function generates samples from a Gaussian copula with two variables,
#' given the specified correlation coefficient between the variables.
#'
#' @param n Integer. The number of samples to generate.
#' @param p Numeric. The correlation coefficient (\eqn{\rho}) between the two variables.
#'          Must be in the range \code{[-1, 1]}.
#' @return A matrix of size \code{n x 2}, where each row represents a sample,
#'         and each column corresponds to one of the two variables. The values
#'         are uniformly distributed in \code{[0, 1]}.
#' @details
#' The function internally constructs a correlation matrix for two variables:
#' \deqn{\rho_matrix = \begin{bmatrix} 1 & p \\ p & 1 \end{bmatrix}}
#' It then calls \code{generate_gaussian_copula_samples} to generate samples.
#'
#' @examples
#' # Example usage:
#' samples <- gaussian_copula_two_vars(n = 1000, p = 0.7)
#' head(samples)
#'
#' @seealso \code{\link{generate_gaussian_copula_samples}}
#' @export
gaussian_copula_two_vars <- function(n, p) {
  # Validate inputs
  if (p < -1 || p > 1) {
    stop("The correlation coefficient 'p' must be between -1 and 1.")
  }

  # Step 1: Define the correlation matrix for two variables
  rho_matrix <- matrix(c(1, p, p, 1), nrow = 2, ncol = 2)

  # Step 2: Generate Gaussian copula samples
  smpls <- generate_gaussian_copula_samples(n, 2, rho_matrix)

  return(smpls)
}
