#' @title Calculate Total Variation (TV) Distance Empirically
#'
#' @description This function calculates the Total Variation (TV) distance between the empirical cumulative distribution
#' functions (ECDFs) of two datasets: original data and generated data. The TV distance is defined as half the sum of
#' the absolute differences between the two CDFs at each point in the domain.
#'
#' @param original_data A numeric vector of the original data.
#' @param generated_data A numeric vector of the generated data.
#' @return A numeric value representing the Total Variation distance between the empirical CDFs of the original and generated data.
#'
#' @examples
#' # Test Case 1: Data from similar distributions
#' original_data <- rnorm(1000, mean = 0, sd = 1)  # Normal distribution (mean = 0, sd = 1)
#' generated_data <- rnorm(1000, mean = 0, sd = 1)  # Similar normal distribution
#' tv_distance <- calculate_tv_distance_empirical(original_data, generated_data)
#' print(tv_distance)  # Expected to be close to 0, as both datasets are similar
#'
#' # Test Case 2: Data from different distributions
#' original_data <- rnorm(1000, mean = 0, sd = 1)  # Normal distribution (mean = 0, sd = 1)
#' generated_data <- rnorm(1000, mean = 5, sd = 2)  # Different normal distribution
#' tv_distance <- calculate_tv_distance_empirical(original_data, generated_data)
#' print(tv_distance)  # Expected to be larger, as the datasets are quite different
#'
#' @export
calculate_tv_distance_empirical <- function(original_data, generated_data) {
  # Create empirical CDFs for original and generated data
  original_ecdf <- ecdf(original_data)
  generated_ecdf <- ecdf(generated_data)

  # Define the grid over which to calculate the distance
  x_values <- sort(unique(c(original_data, generated_data)))

  # Calculate the TV distance as the maximum of absolute differences between ECDFs
  tv_distance <- max(abs(original_ecdf(x_values) - generated_ecdf(x_values)))

  return(tv_distance)
}
