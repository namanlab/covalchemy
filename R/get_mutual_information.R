#' Calculate Mutual Information
#'
#' This function calculates the mutual information between two variables based on their joint distribution.
#' Mutual information measures the amount of information obtained about one variable through the other.
#'
#' @param table A numeric matrix or table. A contingency table or frequency table of two variables.
#' @return A numeric value representing the mutual information between the two variables.
#' @details
#' The mutual information is calculated using the formula:
#' \[ I(X, Y) = H(X) + H(Y) - H(X, Y) \]
#' where:
#' - \( H(X) \) is the entropy of variable X,
#' - \( H(Y) \) is the entropy of variable Y, and
#' - \( H(X, Y) \) is the joint entropy of X and Y.
#'
#' @examples
#' # Example usage with a simple contingency table:
#' pair_table <- table(c(1, 2, 2, 3), c(1, 1, 2, 2))
#' get_mutual_information(pair_table)
#'
#' @export
get_mutual_information <- function(table) {
  # Step 1: Calculate the row sums and column sums
  rowsums <- rowSums(table)
  colsums <- colSums(table)

  # Step 2: Calculate the entropy of each individual variable (X and Y)
  ent_x <- sum(-rowsums / sum(rowsums) * log2(rowsums / sum(rowsums)), na.rm = TRUE)
  ent_y <- sum(-colsums / sum(colsums) * log2(colsums / sum(colsums)), na.rm = TRUE)

  # Step 3: Calculate the mutual information using the formula:
  # I(X, Y) = H(X) + H(Y) - H(X, Y)
  return(ent_x + ent_y - entropy_pair(table))
}
