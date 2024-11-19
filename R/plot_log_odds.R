#' Plot Log-Odds Before and After Transformation
#'
#' This function calculates the log-odds ratio before and after applying a transformation
#' to multiple matrices, and generates a bar plot comparing the log-odds values.
#' The log-odds are calculated using a specified function (default is `log_odds_dc`).
#'
#' @param matrices A list of matrices for which the log-odds are calculated before the transformation.
#' @param new_matrices A list of matrices for which the log-odds are calculated after the transformation.
#' @param names_matrices A vector of names corresponding to the matrices in `matrices` and `new_matrices`.
#' @param log_odds_general A function used to calculate the log-odds (default is `log_odds_dc`).
#'
#' @return A bar plot showing the log-odds before and after the transformation for each matrix and overall.
#'
#' @examples
#' # Example matrices and names
#' matrices <- list(matrix(c(1, 2, 3, 4), nrow = 2), matrix(c(2, 3, 4, 5), nrow = 2))
#' new_matrices <- list(matrix(c(5, 6, 7, 8), nrow = 2), matrix(c(4, 5, 6, 7), nrow = 2))
#' names_matrices <- c("Matrix1", "Matrix2")
#' plot_log_odds(matrices, new_matrices, names_matrices)
#'
#' @import ggplot2 dplyr
#' @export
plot_log_odds <- function(matrices, new_matrices, names_matrices, log_odds_general = log_odds_dc) {

  type <- NULL
  matrix_name <- NULL
  log_odds <- NULL

  # Calculate log-odds for matrices before and after the transformation
  before_log_odds <- sapply(matrices, log_odds_general)
  after_log_odds <- sapply(new_matrices, log_odds_general)

  # Append the name "Overall" to the list of matrix names
  names_matrices <- c(names_matrices, "Overall")

  # Calculate log-odds for the overall matrix (sum of all matrices)
  overall_before <- log_odds_general(Reduce("+", matrices))
  overall_after <- log_odds_general(Reduce("+", new_matrices))

  # Append the overall log-odds to the before and after log-odds vectors
  before_log_odds <- c(before_log_odds, overall_before)
  after_log_odds <- c(after_log_odds, overall_after)

  # Create a data frame for plotting
  data <- data.frame(
    matrix_name = c(names_matrices, names_matrices),
    log_odds = c(before_log_odds, after_log_odds),
    type = c(rep("Before", length(names_matrices)), rep("After", length(names_matrices)))
  )

  # Generate the bar plot using ggplot2
  plot_gen <- data %>%
    mutate(type = factor(type, levels = c("Before", "After"), ordered = TRUE)) %>%
    ggplot(aes(x = matrix_name, y = log_odds, fill = type)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    xlab("Matrix") +
    ylab("Log Odds") +
    ggtitle("Log Odds Before and After Transformation")

  # Print the plot
  print(plot_gen)
}
