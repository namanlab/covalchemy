#' Generate t-Copula Samples
#'
#' This function generates samples from a t-copula given a specified correlation
#' matrix and degrees of freedom. The samples are uniformly distributed in \code{[0, 1]}.
#' across dimensions.
#'
#' @param n Integer. The number of samples to generate.
#' @param d Integer. The dimensionality of the copula.
#' @param rho_matrix A \code{d x d} positive-definite correlation matrix.
#' @param df Numeric. The degrees of freedom of the t-distribution. Must be positive.
#' @return A matrix of size \code{n x d}, where each row represents a sample
#'         and each column corresponds to a dimension. The values are uniformly
#'         distributed in \code{[0, 1]}.
#' @details
#' The function works as follows:
#' \enumerate{
#'   \item Generates multivariate t-distributed samples using the specified
#'         correlation matrix and degrees of freedom.
#'   \item Transforms the samples to the uniform distribution \code{[0, 1]} using the
#'         cumulative distribution function (CDF) of the t-distribution.
#' }
#' @examples
#' # Example usage:
#' library(mvtnorm)  # Load package for `rmvt`
#' rho_matrix <- diag(3)  # 3x3 identity matrix as correlation matrix
#' samples <- generate_t_copula_samples(n = 1000, d = 3, rho_matrix = rho_matrix, df = 5)
#' head(samples)
#'
#' @importFrom mvtnorm rmvt
#' @importFrom stats pt
#' @importFrom stats approxfun coef cor ecdf lm predict quantile rnorm runif
#' @export
generate_t_copula_samples <- function(n, d, rho_matrix, df) {
  # Validate inputs
  if (df <= 0) {
    stop("The degrees of freedom 'df' must be positive.")
  }

  # Step 1: Generate multivariate t-distributed samples
  mvn_samples <- mvtnorm::rmvt(n = n, sigma = rho_matrix, df = df)

  # Step 2: Transform samples to uniform using the CDF of the t-distribution
  uniform_samples <- stats::pt(mvn_samples, df = df)

  return(uniform_samples)
}
