#' Log-Odds Calculation for Concordant and Discordant Pairs
#'
#' This function calculates the log-odds ratio for concordant and discordant pairs
#' based on the contingency table provided. The log-odds ratio is defined as the
#' natural logarithm of the ratio of concordant to discordant pairs.
#'
#' @param tab A contingency table (matrix or data frame) containing counts of pairs
#'            for each combination of outcomes.
#'
#' @return The log-odds ratio calculated as the natural logarithm of the ratio of
#'         concordant pairs to discordant pairs. If discordant pairs are zero,
#'         it returns `Inf` to avoid division by zero.
#'
#' @examples
#' # Example contingency table
#' tab <- matrix(c(10, 5, 7, 8), nrow = 2)
#' log_odds_dc(tab)
#'
#' @importFrom DescTools ConDisPairs
#' @export
log_odds_dc <- function(tab) {

  # Retrieve the concordant and discordant counts from the contingency table
  res <- ConDisPairs(tab)[c("C", "D")]

  # Add 1 to concordant and discordant counts to avoid zero counts
  concordant <- res$C + 1
  discordant <- res$D + 1

  # Check if there are no discordant pairs (discordant == 0)
  if (discordant == 0) {
    return(Inf)  # Return infinity if discordant pairs are zero
  }

  # Calculate and return the log-odds ratio
  return(log(concordant / discordant))
}
