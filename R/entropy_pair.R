#' Calculate Entropy of a Pair
#'
#' This function calculates the entropy of a pair of variables (or a pairwise contingency table)
#' based on the probability distribution of their joint occurrences.
#'
#' @param table A numeric vector or matrix. A contingency table or frequency table of a pair of variables.
#' @return A numeric value representing the entropy of the pair.
#' @details
#' The entropy is calculated using the formula:
#' \[ H(X, Y) = - \sum p(x, y) \cdot \log_2(p(x, y)) \]
#' where \( p(x, y) \) is the probability of observing the pair \( (x, y) \) from the table.
#'
#' @examples
#' # Example usage with a simple contingency table:
#' pair_table <- table(c(1, 2, 2, 3), c(1, 1, 2, 2))
#' entropy_pair(pair_table)
#'
#' @export
entropy_pair <- function(table) {
  # Step 1: Calculate the total number of observations
  total <- sum(table)

  # Step 2: Calculate the probability of each combination in the table
  probs <- table / total

  # Step 3: Calculate the entropy using the probabilities
  # Sum of -p * log2(p) for all pairs (handling NA values with na.rm = TRUE)
  return(sum(-probs * log2(probs), na.rm = TRUE))
}
