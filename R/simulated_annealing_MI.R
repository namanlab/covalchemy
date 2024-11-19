#' Simulated Annealing Algorithm with Target Entropy Stopping Condition
#'
#' This function performs simulated annealing to optimize a given objective (entropy, mutual information, etc.)
#' using a given table modification function. The optimization stops once the target entropy is reached or
#' after a maximum number of iterations.
#'
#' @param initial_table A contingency table to start the optimization.
#' @param obj An objective function that calculates the value to be optimized (e.g., entropy, mutual information).
#' @param gen_fn A function that generates a new table based on the current table.
#' @param target The target value for the objective function (e.g., target entropy).
#' @param max_n The maximum number of iterations to run the algorithm (default is 5000).
#' @param temp The initial temperature for the simulated annealing process (default is 10).
#' @param maxim Logical: Should the algorithm maximize (TRUE) or minimize (FALSE) the objective function (default is TRUE).
#' @param readj Logical: If TRUE, the algorithm is in a readjusting state (default is FALSE).
#'
#' @return A list containing:
#'   - `best`: The best table found during the optimization process.
#'   - `best_eval`: The best evaluation value (objective function value).
#'   - `n`: The number of iterations completed.
#'   - `mutual_info_history`: A data frame with the history of mutual information values (or objective function values)
#'     during each iteration.
#'
#' @examples
#' # Example of using simulated annealing for entropy maximization:
#' initial_table <- matrix(c(5, 3, 4, 2), nrow = 2, ncol = 2)
#' obj <- entropy_pair  # Example entropy function
#' gen_fn <- gen_number_max  # Example generation function
#' target <- 0.5
#' result <- simulated_annealing_MI(initial_table, obj, gen_fn, target)
#'
#' @export
simulated_annealing_MI <- function(initial_table, obj, gen_fn, target,
                         max_n = 5000, temp = 10, maxim = TRUE, readj = FALSE) {

  # Step 1: Initialize best and current solutions
  best <- initial_table
  best_eval <- obj(best)  # Evaluate the objective for the best solution
  curr <- best
  curr_eval <- best_eval  # Evaluate the objective for the current solution
  n <- 0  # Iteration counter

  # Initialize the history of mutual information (or objective function values)
  mutual_info_history <- data.frame(iteration = 0, mutual_info = curr_eval,
                                    type = ifelse(readj, "Readjusting", ifelse(maxim, "Maximizing", "Minimizing")))

  # Step 2: Start the annealing loop
  while (n < max_n) {

    # Check if the target entropy (or objective) is reached
    if ((maxim && curr_eval >= target) || (!maxim && curr_eval <= target)) {
      print(paste("Target entropy reached:", curr_eval))
      break
    }

    # Generate a candidate solution
    cand <- gen_fn(curr)
    cand_eval <- obj(cand)  # Evaluate the candidate solution

    # Step 3: Determine whether to accept the candidate based on the objective function
    if (maxim) {
      # Maximize the objective
      if (cand_eval > best_eval) {
        best <- cand
        best_eval <- cand_eval  # Update best if candidate is better
      }
      diff <- cand_eval - curr_eval
    } else {
      # Minimize the objective
      if (cand_eval < best_eval) {
        best <- cand
        best_eval <- cand_eval  # Update best if candidate is better
      }
      diff <- -cand_eval + curr_eval
    }

    # Step 4: Simulated annealing temperature schedule
    temp <- 0.9 * temp  # Decrease temperature over time
    metropolis <- exp(diff / temp)  # Metropolis acceptance probability

    # Step 5: Decide whether to accept the candidate solution
    if (diff > 0 || runif(1) < metropolis) {
      curr <- cand
      curr_eval <- cand_eval  # Accept the candidate
    }

    # Step 6: Record mutual information history
    n <- n + 1  # Increment the iteration counter
    mutual_info_history <- rbind(mutual_info_history, data.frame(iteration = n, mutual_info = curr_eval,
                                                                 type = ifelse(readj, "Readjusting",
                                                                               ifelse(maxim, "Maximizing", "Minimizing"))))
  }

  # Return results: best table, best evaluation, iteration count, and history
  return(list(best = best, best_eval = best_eval, n = n, mutual_info_history = mutual_info_history))
}
