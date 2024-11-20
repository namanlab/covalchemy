#' @title Get Optimal Grid Assignment
#'
#' @description This function computes an optimal grid assignment between two variables `x` and `y` based on a third variable `z`.
#' It uses the quantiles of `x` and `y` to segment the data based on the distribution of `z`. Then, it computes the cost of
#' assigning points from `x` to `y` by calculating the counts of `y` values within quantile ranges of `x` and `y`, and then solves
#' the assignment problem using the Hungarian algorithm.
#'
#' @param x A numeric vector of values representing the first variable.
#' @param y A numeric vector of values representing the second variable.
#' @param z A numeric vector representing the distribution of the third variable, used to define quantiles for `x` and `y`.
#' @return A numeric vector of optimal indices that represents the optimal assignment of `x` values to `y` values.
#'
#' @examples
#' # Test Case 1: Simple uniform data
#' x <- rnorm(1000)
#' y <- rnorm(1000)
#' z <- sample(1:5, 1000, replace = TRUE)
#' optimal_assignment <- get_optimal_grid(x, y, z)
#'
#' # Test Case 2: Data with a skewed distribution
#' x <- rexp(1000, rate = 1)
#' y <- rpois(1000, lambda = 2)
#' z <- sample(1:3, 1000, replace = TRUE)
#' optimal_assignment <- get_optimal_grid(x, y, z)
#'
#' @importFrom clue solve_LSAP
#' @export
get_optimal_grid <- function(x, y, z) {
  # Calculate quantiles for x and y based on z distribution
  z_levels <- sort(unique(z))  # Get unique levels of z
  x_quantiles <- quantile(x, probs = cumsum(table(z) / length(z)))  # Quantiles for x based on z distribution
  y_quantiles <- quantile(y, probs = cumsum(table(z) / length(z)))  # Quantiles for y based on z distribution

  k = length(z_levels)  # Number of unique levels in z

  # Initialize result matrix for storing counts of assignments
  res = matrix(rep(0, k^2), nrow = k)  # Create a k x k matrix initialized with 0

  # Iterate over unique levels of z to fill the matrix
  for (i in seq_along(z_levels)) {
    # Get indices of x within the quantile range for the current z level
    x_range <- if (i == 1) {
      which(x <= x_quantiles[i])  # For the first level, select values <= quantile
    } else {
      which(x > x_quantiles[i - 1] & x <= x_quantiles[i])  # Otherwise, select values within quantile range
    }

    for (j in seq_along(z_levels)) {
      # Filter y values corresponding to the selected x values
      temp_y = y[x_range]

      # Get indices of y within the quantile range for the current z level
      y_range <- if (j == 1) {
        which(temp_y <= y_quantiles[j])  # For the first level, select values <= quantile
      } else {
        which(temp_y > y_quantiles[j - 1] & temp_y <= y_quantiles[j])  # Otherwise, select values within quantile range
      }

      # Store the count of y values in the quantile range in the result matrix
      res[i, j] = length(y_range)
    }
  }

  # Convert the result matrix to a cost matrix by subtracting from the maximum count
  M <- max(res)  # Get the maximum value in the result matrix
  cost_matrix <- M - res  # Cost matrix: max count minus the counts in res

  # Solve the assignment problem using the Hungarian algorithm (Linear Sum Assignment Problem)
  assignment <- solve_LSAP(cost_matrix)

  # Return the row and column indices of the selected cells as a vector
  return(as.vector(assignment))
}
