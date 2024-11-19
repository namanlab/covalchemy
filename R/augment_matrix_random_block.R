#' Augment Matrix with Random 2x2 Block Adjustment
#'
#' This function selects a random 2x2 block of values in the input matrix `table` and modifies them
#' based on the specified delta. It checks certain conditions before applying the modifications to
#' the selected block. The process repeats until a valid block is found or a maximum of 100 iterations
#' is reached.
#'
#' @param table A matrix (numeric) to which the random block adjustment will be applied.
#' @param delta A numeric value that determines the magnitude of the adjustment.
#'              If positive, values are subtracted from the block; if negative, values are added.
#'
#' @return A matrix (numeric) with the adjusted 2x2 block.
#'
#' @examples
#' table <- matrix(1:9, 3, 3)
#' augment_matrix_random_block(table, 1)
#'
#' @export
augment_matrix_random_block <- function(table, delta) {
  nrows <- nrow(table)  # Number of rows in the matrix
  ncols <- ncol(table)  # Number of columns in the matrix
  iters <- 0             # Initialize iteration counter

  repeat {
    iters <- iters + 1  # Increment iteration counter

    # Randomly select two rows and two columns
    RR <- sort(sample(1:nrows, 2))  # Row indices for the block
    CC <- sort(sample(1:ncols, 2))  # Column indices for the block

    # Check conditions before applying delta
    if (delta < 0 && table[RR[1], CC[1]] > 0 && table[RR[2], CC[2]] > 0) {
      break  # If delta is negative and conditions are met, break out of loop
    }
    if (delta > 0 && table[RR[1], CC[2]] > 0 && table[RR[2], CC[1]] > 0) {
      break  # If delta is positive and conditions are met, break out of loop
    }

    # If iteration count exceeds 100, set delta to 0 and break
    if (iters > 100) {
      delta <- 0
      break
    }
  }

  # Modify the selected block of the matrix based on delta
  table[RR[1], CC[1]] <- table[RR[1], CC[1]] + delta  # Add delta to top-left
  table[RR[2], CC[2]] <- table[RR[2], CC[2]] + delta  # Add delta to bottom-right
  table[RR[1], CC[2]] <- table[RR[1], CC[2]] - delta  # Subtract delta to top-right
  table[RR[2], CC[1]] <- table[RR[2], CC[1]] - delta  # Subtract delta to bottom-left

  return(table)  # Return the modified table
}
