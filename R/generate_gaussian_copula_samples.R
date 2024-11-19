#' Generate Gaussian Copula Samples
#'
#' This function generates samples from a Gaussian copula given a specified
#' correlation matrix. The samples are uniformly distributed in \code{[0, 1]} across
#' dimensions.
#'
#' @param n Integer. The number of samples to generate.
#' @param d Integer. The dimensionality of the copula.
#' @param rho_matrix A \code{d x d} positive-definite correlation matrix.
#' @return A matrix of size \code{n x d}, where each row represents a sample
#'         and each column corresponds to a dimension. The values are uniformly
#'         distributed in \code{[0, 1]}.
#' @details
#' The function works as follows:
#' \enumerate{
#'   \item Generates multivariate normal samples with the given correlation matrix.
#'   \item Transforms the samples to the uniform distribution \code{[0, 1]} using the
#'         cumulative distribution function (CDF) of the standard normal.
#' }
#' @examples
#' # Example usage:
#' library(MASS)  # Load package for `mvrnorm`
#' rho_matrix <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)  # 2x2 correlation matrix
#' samples <- generate_gaussian_copula_samples(n = 1000, d = 2, rho_matrix = rho_matrix)
#' head(samples)
#'
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats pnorm
#' @importFrom stats approxfun coef cor ecdf lm predict quantile rnorm runif
#' @export
generate_gaussian_copula_samples <- function(n, d, rho_matrix) {
  # Step 1: Generate multivariate normal samples
  mean_vector <- rep(0, d)  # Mean vector for multivariate normal
  mvn_samples <- mvtnorm::rmvnorm(n = n, mean = mean_vector, sigma = rho_matrix)

  # Step 2: Transform samples to uniform using the CDF of the standard normal
  uniform_samples <- stats::pnorm(mvn_samples)

  return(uniform_samples)
}
