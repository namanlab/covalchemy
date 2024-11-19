#' Generate Samples with Target Kendall's Tau Correlation Using a Copula Approach
#'
#' This function generates two variables with a specified target Kendall's tau correlation
#' using copula-based methods. The user can specify the type of copula (Gaussian or t),
#' the type of inverse CDF method to apply to the variables, and the degree of the polynomial
#' interpolation if applicable.
#'
#' @param x1 A numeric vector. The first dataset used for generating inverse CDFs.
#' @param x2 A numeric vector. The second dataset used for generating inverse CDFs.
#' @param target_corr_kendall A numeric value. The desired target Kendall's tau correlation
#'                            between the two generated variables.
#' @param copula_type A string. The type of copula to use, either "gaussian" or "t" (default is "gaussian").
#' @param inv_cdf_type A string. The type of inverse CDF method to use. Options include:
#'                     "quantile_1", "quantile_4", "quantile_7", "quantile_8", "linear",
#'                     "akima", "poly" (default is "quantile_7").
#' @param degree An integer. The degree of the polynomial interpolation (default is 10).
#' @return A list containing two components: \code{x1} and \code{x2}, which are the modified
#'         versions of the input datasets \code{x1} and \code{x2} with the desired target
#'         Kendall's tau correlation.
#' @details
#' This function works by:
#' 1. Generating two variables using the specified copula type (Gaussian or t) with the target
#'    Kendall's tau correlation.
#' 2. Applying the chosen inverse CDF transformation to the generated copula samples.
#' 3. Returning the modified variables that have the target correlation.
#'
#' @examples
#' # Example usage:
#' x1 <- ChickWeight$weight
#' x2 <- ChickWeight$Time
#' cor(x1, x2, method = "kendall")  # Calculate original Kendall's tau correlation
#' res <- get_target_corr(x1, x2, target_corr_kendall = 0,
#'                        copula_type = "gaussian", inv_cdf_type = "poly")
#' cor(res$x1, res$x2, method = "kendall")  # Calculate modified Kendall's tau correlation
#'
#' @seealso \code{\link{gaussian_copula_two_vars}}, \code{\link{t_copula_two_vars}},
#'          \code{\link{genCDFInv_quantile}}, \code{\link{genCDFInv_linear}},
#'          \code{\link{genCDFInv_akima}}, \code{\link{genCDFInv_poly}}
#' @export
get_target_corr <- function(x1, x2, target_corr_kendall,
                            copula_type = "gaussian", inv_cdf_type = "quantile_7",
                            degree = 10) {
  # Step 1: Generate copula samples with target correlation
  n <- length(x1)

  if (copula_type == "gaussian") {
    target_rho <- sin(target_corr_kendall * pi / 2)  # Using Greiner's equality for Gaussian copula
    res <- gaussian_copula_two_vars(n, target_rho)
  } else if (copula_type == "t") {
    target_rho <- sin(target_corr_kendall * pi / 2)  # Using Greiner's equality for t-copula
    res <- t_copula_two_vars(n, target_rho)
  } else {
    stop("Unsupported copula type")
  }

  # Step 2: Apply chosen inverse CDF transformation to both variables
  inv_cdf_d1 <- switch(inv_cdf_type,
                       "quantile_1" = genCDFInv_quantile(x1, 1),   # Stepwise inverse CDF
                       "quantile_4" = genCDFInv_quantile(x1, 4),   # Linear interpolation
                       "quantile_7" = genCDFInv_quantile(x1, 7),   # Default
                       "quantile_8" = genCDFInv_quantile(x1, 8),   # Median-unbiased
                       "linear" = genCDFInv_linear(x1),             # Linear interpolation
                       "akima" = genCDFInv_akima(x1),               # Akima spline interpolation
                       "poly" = genCDFInv_poly(x1, degree = degree))  # Polynomial regression

  inv_cdf_d2 <- switch(inv_cdf_type,
                       "quantile_1" = genCDFInv_quantile(x2, 1),   # Stepwise inverse CDF
                       "quantile_4" = genCDFInv_quantile(x2, 4),   # Linear interpolation
                       "quantile_7" = genCDFInv_quantile(x2, 7),   # Default
                       "quantile_8" = genCDFInv_quantile(x2, 8),   # Median-unbiased
                       "linear" = genCDFInv_linear(x2),             # Linear interpolation
                       "akima" = genCDFInv_akima(x2),               # Akima spline interpolation
                       "poly" = genCDFInv_poly(x2, degree = degree))  # Polynomial regression

  # Step 3: Vectorize the inverse CDF functions
  F1Inv <- Vectorize(inv_cdf_d1)
  F2Inv <- Vectorize(inv_cdf_d2)

  # Step 4: Apply the inverse CDF transformations to the copula samples
  x1mod <- F1Inv(res[, 1])
  x2mod <- F2Inv(res[, 2])

  return(list(x1 = x1mod, x2 = x2mod))
}
