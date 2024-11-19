#' @title Get Target Entropy
#'
#' @description This function adjusts the mutual information between two categorical variables (x and y)
#' by modifying their contingency table using simulated annealing to reach a target entropy.
#' The function first calculates the range of possible entropy values (min and max) and checks if the target entropy
#' lies within that range. If so, it adjusts the mutual information to reach the target entropy, either by increasing or
#' decreasing it, depending on the initial entropy value.
#'
#' @param x A vector of categorical values representing variable x.
#' @param y A vector of categorical values representing variable y.
#' @param target_entropy The target entropy value to be reached.
#' @param max_n Maximum number of iterations for the optimization process (default is 10,000).
#' @param epsilon The tolerance value for determining if the target entropy has been reached (default is 0.001).
#' @return A list containing:
#'   - final_df: A dataframe with the adjusted contingency table.
#'   - final_table: The final contingency table after adjustments.
#'   - history: The history of the optimization process.
#'   - max_mut: The maximum mutual information found.
#'   - min_mut: The minimum mutual information found.
#'
#' @examples
#' set.seed(33)
#' df <- data.frame(
#'   x = sample(paste("Categ", 1:4), 10000, replace = TRUE),
#'   y = sample(paste("Categ", 10:4), 10000, replace = TRUE)
#' )
#' target_entropy <- 1  # Set your target entropy here
#' \dontrun{res <- get_target_entropy(df$x, df$y, target_entropy)}
#'
#'
#' @export
get_target_entropy <- function(x, y, target_entropy, max_n = 10000, epsilon = 0.001) {

  # Create contingency table
  table <- table(x, y)

  # Get initial mutual information
  old_mut <- get_mutual_information(table)

  # Step 1: Find range of entropy values
  # Max Entropy (Maximizing mutual information)
  result_max <- simulated_annealing_MI(table, obj = get_mutual_information, gen_fn = gen_number_max,
                             target = Inf, max_n = max_n, maxim = TRUE)
  new_mut_max <- result_max[[2]]

  # Min Entropy (Minimizing mutual information)
  result_min_sk <- sinkhorn_algorithm(table, obj = get_mutual_information, max_iter = max_n)
  new_mut_min <- result_min_sk[[2]]

  # Print the max and min entropy values
  print(paste("Current Entropy:", old_mut))
  print(paste("Min Entropy:", new_mut_min))
  print(paste("Max Entropy:", new_mut_max))

  # Step 2: Check if the target entropy is within the range
  if (target_entropy < new_mut_min || target_entropy > new_mut_max) {
    print("Target entropy is out of range. Please choose a value between the min and max entropy.")
    return(NULL)
  }

  # Step 3: Adjust mutual information to reach the target entropy
  if (target_entropy > old_mut) {
    gen_fn <- gen_number_max
    result <- simulated_annealing_MI(table, obj = get_mutual_information, gen_fn = gen_fn,
                           target = target_entropy, max_n = max_n, maxim = TRUE)
    final_hist = result[[4]]
    final_table <- result[[1]]
    final_mut <- result[[2]]
    if (result[[2]] - target_entropy > epsilon) {
      print("Target exceeded. Re-adjusting with gen_number_1 to decrease.")
      result_sub <- simulated_annealing_MI(result[[1]], obj = get_mutual_information, gen_fn = gen_number_1,
                                 target = target_entropy, max_n = max_n,
                                 maxim = FALSE, readj = T)
      result_sub[[4]]$iteration = result_sub[[4]]$iteration + max(final_hist$iteration)
      final_hist = rbind(result_sub[[4]], final_hist)
      final_table <- result_sub[[1]]
      final_mut <- result_sub[[2]]
    }

  } else {
    gen_fn <- gen_number_min
    result <- simulated_annealing_MI(table, obj = get_mutual_information, gen_fn = gen_fn,
                           target = target_entropy, max_n = max_n, maxim = FALSE)
    final_hist = result[[4]]
    final_table <- result[[1]]
    final_mut <- result[[2]]
    if (target_entropy - result[[2]] > epsilon ) {
      print("Target crossed. Re-adjusting with gen_number_1 to increase.")
      result_sub <- simulated_annealing_MI(result[[1]], obj = get_mutual_information, gen_fn = gen_number_1,
                                 target = target_entropy, max_n = max_n,
                                 maxim = TRUE, readj = T)
      result_sub[[4]]$iteration = result_sub[[4]]$iteration + max(final_hist$iteration)
      final_hist = rbind(result_sub[[4]], final_hist)
      final_table <- result_sub[[1]]
      final_mut <- result_sub[[2]]
    }
  }

  # Print the final mutual information value
  print(paste("Final Mutual Information:", final_mut))

  # Return results
  contingency_df <- as.data.frame(as.table(final_table))
  expanded_df <- contingency_df[rep(1:nrow(contingency_df), contingency_df$Freq), c("x", "y")]

  return(list(final_df = expanded_df, final_table = final_table, history = final_hist,
              max_mut = new_mut_max, min_mut = new_mut_min))
}
