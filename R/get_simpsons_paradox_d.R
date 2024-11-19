#' Introduce Simpson's Paradox in Discrete Data
#'
#' This function modifies contingency tables associated with different levels of a categorical variable
#' to create or highlight Simpson's Paradox using simulated annealing. The paradox occurs when aggregated
#' data trends differ from subgroup trends.
#'
#' @param x A vector of categorical values for the first variable.
#' @param y A vector of categorical values for the second variable.
#' @param z A vector indicating levels of a third variable that segments the data.
#' @param manual_vec A numeric vector specifying target log-odds trends for each level of `z`.
#' @param target_overall A numeric value representing the target log-odds for the aggregated data.
#' @param margin A numeric value for allowed deviation in log-odds within each subgroup.
#' @param margin_overall A numeric value for allowed deviation in aggregated log-odds.
#' @param max_n An integer specifying the maximum number of iterations for the annealing process.
#' @param temp A numeric value for the initial temperature in the annealing process.
#' @param log_odds_general A function to compute the log-odds for a given contingency table (default: `log_odds_dc`).
#'
#' @return A list containing:
#' - `final_df`: A data frame representing the modified dataset.
#' - `final_table`: A list of modified contingency tables.
#' - `history`: A data frame tracking the overall log-odds over iterations.
#'
#' @details
#' This function works by iteratively modifying individual matrices (contingency tables) corresponding
#' to levels of `z` while respecting log-odds constraints. The overall log-odds of the aggregated table
#' are also adjusted to achieve the specified `target_overall`. Simulated annealing ensures that the
#' modifications balance between achieving the targets and avoiding overfitting.
#'
#' @examples
#' # Example with predefined contingency tables
#' set.seed(42)
#' matrices <- list(
#'   ta = matrix(c(512, 89, 313, 19), ncol = 2, byrow = TRUE),
#'   tb = matrix(c(353, 17, 207, 8), ncol = 2, byrow = TRUE),
#'   tc = matrix(c(120, 202, 205, 391), ncol = 2, byrow = TRUE)
#' )
#' df_list <- lapply(seq_along(matrices), function(i) {
#'   mat <- matrices[[i]]
#'   z_level <- names(matrices)[i]
#'   df <- as.data.frame(as.table(mat))
#'   colnames(df) <- c("x", "y", "Freq")
#'   df$z <- z_level
#'   return(df)
#' })
#' final_df <- do.call(rbind, df_list)
#' expanded_df <- final_df[rep(1:nrow(final_df), final_df$Freq), c("x", "y", "z")]
#' result <- get_simpsons_paradox_d(
#'   expanded_df$x, expanded_df$y, expanded_df$z,
#'   manual_vec = c(-1, -1, -1),
#'   target_overall = +1,
#'   margin = 0.2, margin_overall = 0.2, max_n = 200
#' )
#' table(expanded_df$x) - table(result$final_df$x)
#'
#' @export
get_simpsons_paradox_d <- function(x, y, z,
                                   manual_vec, target_overall,
                                   margin, margin_overall,
                                   max_n = 1000, temp = 10,
                                   log_odds_general = log_odds_dc) {

  # Helper function to sum matrices
  sum_matrices <- function(mats) Reduce("+", mats)

  # Generate matrices (contingency tables) based on levels of z
  levels_z <- unique(z)
  matrices <- lapply(levels_z, function(level) {
    # Subset the data for the current level of z
    subset_data <- data.frame(x = x[z == level], y = y[z == level])
    # Create the contingency table for the current level
    table(subset_data$x, subset_data$y)
  })
  names(matrices) <- levels_z

  curr <- matrices
  curr_eval_sum <- log_odds_general(sum_matrices(curr))
  history <- data.frame(iteration = 0, overall_log_odds = curr_eval_sum)
  for (n in 1:max_n) {
    log_odds_values <- sapply(curr, log_odds_general)
    selection_probs <- sapply(1:length(curr), function(i) {
      cur_log_odds <- log_odds_values[i]
      target <- manual_vec[i]
      if (cur_log_odds < -margin && target < 0) return(-Inf)  # No modification
      if (cur_log_odds > margin && target > 0) return(-Inf)  # No modification
      return(abs(cur_log_odds - target))
    })
    # print(selection_probs)
    if (all(is.infinite(selection_probs))) {selection_probs <- rep(0, length(selection_probs))}
    selection_probs <- softmax(selection_probs)  # Normalize to softmax probabilities
    idx <- sample(1:length(curr), 1, prob = selection_probs)
    cur_log_odds <- log_odds_general(curr[[idx]])
    target <- manual_vec[idx]
    delta <- 1
    if (cur_log_odds - target*margin > 0){delta = -1}
    modified_matrix <- augment_matrix_random_block(curr[[idx]], delta)
    new_matrices <- curr
    new_matrices[[idx]] <- modified_matrix
    new_eval_sum <- log_odds_general(sum_matrices(new_matrices))
    temp <- 0.9 * temp  # Annealing step
    metropolis <- ifelse(delta <= 0, exp((curr_eval_sum - new_eval_sum) / temp),
                         exp((new_eval_sum - curr_eval_sum)/temp))
    # Metropolis criteria for acceptance and constraints on log-odds limits
    if (runif(1) < metropolis) {
      curr <- new_matrices
      curr_eval_sum <- new_eval_sum
    }
    history <- rbind(history, data.frame(iteration = n, overall_log_odds = curr_eval_sum))
  }

  ## NOW ADJUST OVERALL
  curr_eval_sum <- log_odds_general(sum_matrices(curr))
  if (curr_eval_sum < margin_overall*target_overall && target_overall > 0){
    # increase
    for (n in 1:max_n) {
      log_odds_values <- sapply(curr, log_odds_general)
      selection_probs <- sapply(1:length(curr), function(i) {
        cur_log_odds <- log_odds_values[i]
        target <- manual_vec[i]
        if (target > 0){return(1000*exp(-cur_log_odds))}
        else {return(-margin - cur_log_odds)}
      })
      if (all(is.infinite(selection_probs))) {selection_probs <- rep(0, length(selection_probs))}
      selection_probs <- softmax(selection_probs)  # Normalize to softmax probabilities
      idx <- sample(1:length(curr), 1, prob = selection_probs)
      delta <- 1
      modified_matrix <- augment_matrix_random_block(curr[[idx]], delta)
      new_matrices <- curr
      new_matrices[[idx]] <- modified_matrix
      new_eval_sum <- log_odds_general(sum_matrices(new_matrices))
      temp <- 0.9 * temp  # Annealing step
      metropolis <- ifelse(delta <= 0, exp((curr_eval_sum - new_eval_sum) / temp),
                           exp((new_eval_sum - curr_eval_sum)/temp))
      # CHECK IF MODIFIED MATRIX CONSTRAINTS ARE NOT VIOLATED!
      cur_log_odds <- log_odds_general(modified_matrix)
      target <- manual_vec[idx]
      CHECK_F = 0
      if (cur_log_odds < -margin && target < 0){CHECK_F = 1}
      if (cur_log_odds > margin && target > 0){CHECK_F = 1}
      # Metropolis criteria for acceptance and constraints on log-odds limits
      if (runif(1) < metropolis && CHECK_F == 1) {
        curr <- new_matrices
        curr_eval_sum <- new_eval_sum
      }
      history <- rbind(history, data.frame(iteration = n + max_n, overall_log_odds = curr_eval_sum))
    }
  } else if (curr_eval_sum > margin_overall*target_overall && target_overall < 0){
    # decrease
    for (n in 1:max_n) {
      log_odds_values <- sapply(curr, log_odds_general)
      selection_probs <- sapply(1:length(curr), function(i) {
        cur_log_odds <- log_odds_values[i]
        target <- manual_vec[i]
        if (target < 0){return(1000*exp(cur_log_odds))}
        else {return(cur_log_odds - margin)}
      })
      if (all(is.infinite(selection_probs))) {selection_probs <- rep(0, length(selection_probs))}
      selection_probs <- softmax(selection_probs)  # Normalize to softmax probabilities
      idx <- sample(1:length(curr), 1, prob = selection_probs)
      delta <- -1
      modified_matrix <- augment_matrix_random_block(curr[[idx]], delta)
      new_matrices <- curr
      new_matrices[[idx]] <- modified_matrix
      new_eval_sum <- log_odds_general(sum_matrices(new_matrices))
      temp <- 0.9 * temp  # Annealing step
      metropolis <- ifelse(delta <= 0, exp((curr_eval_sum - new_eval_sum) / temp),
                           exp((new_eval_sum - curr_eval_sum)/temp))
      # CHECK IF MODIFIED MATRIX CONSTRAINTS ARE NOT VIOLATED!
      cur_log_odds <- log_odds_general(modified_matrix)
      target <- manual_vec[idx]
      CHECK_F = 0
      if (cur_log_odds < -margin && target < 0){CHECK_F = 1}
      if (cur_log_odds > margin && target > 0){CHECK_F = 1}
      # Metropolis criteria for acceptance and constraints on log-odds limits
      if (runif(1) < metropolis && CHECK_F == 1) {
        curr <- new_matrices
        curr_eval_sum <- new_eval_sum
      }
      history <- rbind(history, data.frame(iteration = n + max_n, overall_log_odds = curr_eval_sum))
    }
  }


  # get back dataframe from list of matrices
  df_list <- lapply(seq_along(curr), function(i) {
    mat <- curr[[i]]
    z_level <- levels_z[i]
    df <- as.data.frame(as.table(mat))
    colnames(df) <- c("x", "y", "Freq")
    df$z <- z_level
    return(df)
  })

  # Combine all the data frames into one
  final_df <- do.call(rbind, df_list)
  expanded_df <- final_df[rep(1:nrow(final_df), final_df$Freq), c("x", "y", "z")]
  plot_log_odds(matrices, curr, names(matrices))

  # Return the final matrices and log odds history
  return(list(final_df = expanded_df, final_table = curr,
              history = history))
}
