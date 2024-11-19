#' Generate a New Number for Minimizing Mutual Information
#'
#' This function modifies a given contingency table by swapping values between two cells
#' in a way that minimizes the mutual information. The function randomly selects two cells
#' from the table and adjusts their values to reduce mutual information, returning the modified table.
#'
#' @param x A contingency table (numeric matrix or table).
#' @return A modified contingency table with minimized mutual information.
#' @details
#' This function performs the following steps:
#' 1. Randomly selects two rows and two columns from the table.
#' 2. Ensures that the selected cells have non-zero values.
#' 3. Adjusts the values of the selected cells in the table to minimize mutual information.
#' 4. Returns the modified table with minimized mutual information.
#'
#' @examples
#' # Example usage with a contingency table:
#' pair_table <- table(c(1, 2, 2, 3), c(1, 1, 2, 2))
#' gen_number_min(pair_table)
#'
#' @export
gen_number_min <- function(x) {
  # Step 1: Extract dimensions of the table
  table <- x
  nrows <- nrow(table)
  ncols <- ncol(table)

  repeat {
    # Step 2: Randomly select two rows and two columns
    RR <- sample(1:nrows, 2)  # Two rows
    CC <- sample(1:ncols, 2)  # Two columns

    # Step 3: Ensure that the selected cells have non-zero values
    if (table[RR[1], CC[1]] > 0 && table[RR[2], CC[2]] > 0) break
  }

  # Step 4: Calculate the sum of values of the selected cells
  S <- (table[RR[1], CC[1]] + table[RR[2], CC[2]] + table[RR[1], CC[2]] + table[RR[2], CC[1]])

  # Step 5: Calculate delta to adjust the values for minimizing mutual information
  delta <- round((table[RR[1], CC[1]] * table[RR[2], CC[2]] - table[RR[1], CC[2]] * table[RR[2], CC[1]]) / S)

  # Step 6: Modify the values of the selected cells
  table[RR[1], CC[1]] <- table[RR[1], CC[1]] - delta
  table[RR[2], CC[2]] <- table[RR[2], CC[2]] - delta
  table[RR[1], CC[2]] <- table[RR[1], CC[2]] + delta
  table[RR[2], CC[1]] <- table[RR[2], CC[1]] + delta

  # Step 7: Return the modified table
  return(table)
}
