#' Sinkhorn Algorithm for Matrix Scaling
#'
#' This function applies the Sinkhorn-Knopp algorithm to adjust the row and column sums of a matrix
#' to match the target sums. The algorithm iteratively scales the rows and columns by updating
#' scaling factors (alpha and beta) until convergence or the maximum number of iterations is reached.
#'
#' @param initial_table A matrix to be adjusted using the Sinkhorn algorithm.
#' @param obj An objective function to evaluate the matrix (e.g., entropy, mutual information).
#' @param max_iter The maximum number of iterations for the algorithm (default is 500).
#' @param tolerance The convergence tolerance. If the change in the objective function is smaller than this,
#'        the algorithm stops (default is 1e-5).
#'
#' @return A list containing:
#'   - `updated_table`: The matrix after Sinkhorn scaling.
#'   - `new_mut`: The objective function value for the scaled matrix.
#'   - `iter`: The number of iterations performed.
#'   - `mutual_info_history`: A data frame with the history of objective function values during each iteration.
#'
#' @examples
#' initial_table <- matrix(c(5, 3, 4, 2), nrow = 2, ncol = 2)
#' obj <- entropy_pair  # Example entropy function
#' result <- sinkhorn_algorithm(initial_table, obj)
#'
#' @export
sinkhorn_algorithm <- function(initial_table, obj, max_iter = 500, tolerance = 1e-5) {

  # Step 1: Initialize variables
  S <- sum(initial_table)  # Total sum of the matrix
  S_r <- rowSums(initial_table)  # Row sums
  S_c <- colSums(initial_table)  # Column sums
  n <- nrow(initial_table)  # Number of rows
  m <- ncol(initial_table)  # Number of columns

  # Initialize scaling factors (alpha for columns, beta for rows)
  alpha <- rep(1, m)  # Column scaling factors
  beta <- rep(1, n)   # Row scaling factors

  # Evaluate the initial matrix with the given objective function
  curr_eval <- obj(initial_table)

  # Initialize history of objective values
  mutual_info_history <- data.frame(iteration = 0, mutual_info = curr_eval,
                                    type = "Minimizing")

  # Step 2: Sinkhorn scaling loop
  for (iter in 1:max_iter) {

    # Create a matrix of ones with the same dimensions as the initial matrix
    ones_mat <- matrix(1, nrow = m, ncol = n)

    # Step 3: Update alpha (column scaling)
    alpha <- S_c / (ones_mat %*% beta)

    # Step 4: Update beta (row scaling)
    beta <- S_r / (t(ones_mat) %*% alpha)

    # Step 5: Construct the updated matrix using the scaling factors
    updated_table <- t(diag(as.vector(alpha)) %*% ones_mat %*% diag(as.vector(beta)))

    # Step 6: Evaluate the updated table with the objective function
    curr_eval <- obj(updated_table)

    # Step 7: Record the current evaluation in the history
    mutual_info_history <- rbind(mutual_info_history, data.frame(iteration = iter,
                                                                 mutual_info = curr_eval,
                                                                 type = "Minimizing"))

    # Step 8: Check for convergence based on the tolerance (if needed)
    if (iter > 1 && abs(curr_eval - mutual_info_history$mutual_info[iter - 1]) < tolerance) {
      print(paste("Converged at iteration", iter, "with objective value", curr_eval))
      break
    }
  }

  # Final evaluation of the updated table
  new_mut <- obj(updated_table)

  # Return results: the updated table, final objective value, number of iterations, and history
  return(list(updated_table = updated_table, new_mut = new_mut, iter = iter,
              mutual_info_history = mutual_info_history))
}
