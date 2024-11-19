#' Softmax Function with Special Handling for Infinite Values
#'
#' This function computes the softmax of a vector `x`, with special handling for infinite values.
#' The softmax function transforms input values into a probability distribution by exponentiating
#' each value, then normalizing by the sum of all exponentiated values.
#' The function ensures numerical stability, particularly when dealing with very large or very small values,
#' and handles cases where the values are infinite (`Inf`).
#'
#' @param x A numeric vector for which the softmax function will be calculated.
#'
#' @return A numeric vector of the same length as `x`, where the values represent probabilities summing to 1.
#'
#' @examples
#' softmax(c(10, 5, 2))
#' softmax(c(Inf, -Inf, 0))
#' softmax(c(-Inf, -Inf, -Inf))
#'
#' @export
softmax <- function(x) {
  # Check for any positive infinity values in `x` and handle them
  # If any value is Inf (positive), return 1 for that position and 0 for others
  if (any(is.infinite(x) & x > 0)) {
    probs <- rep(0, length(x))   # Initialize all probabilities to 0
    probs[which(x == Inf)] <- 1  # Set the probability of Inf values to 1
    return(probs)
  }

  # Handle the case where all values in `x` are negative infinity (-Inf)
  if (all(x == -Inf)) {
    return(rep(0, length(x)))  # All values will have probability 0
  }

  # Standard softmax calculation with numerical stability
  # Adjust values to avoid overflow when computing exp(x) for large values
  exp_x <- exp(x - max(x, na.rm = TRUE))  # Subtract max for stability
  res <- exp_x / sum(exp_x, na.rm = TRUE)  # Normalize to get probabilities

  # Return normalized softmax probabilities
  return(res / sum(res))  # Ensures the sum of probabilities equals 1
}
