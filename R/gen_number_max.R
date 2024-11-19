#' Generate a New Number for Maximizing Mutual Information
#'
#' This function modifies a given contingency table by swapping values between two cells
#' to maximize the mutual information. The function randomly selects two cells from the table
#' and adjusts their values in a way that increases mutual information. The function then
#' returns the modified table with the highest mutual information.
#'
#' @param x A contingency table (numeric matrix or table).
#' @return A modified contingency table with maximized mutual information.
#' @details
#' This function performs the following steps:
#' 1. Randomly selects two rows and two columns from the table.
#' 2. Checks if the selected cells have non-zero values.
#' 3. Adjusts the values of the selected cells in two different modified tables, `table1` and `table2`.
#' 4. Calculates the mutual information for both modified tables.
#' 5. Returns the table with the higher mutual information.
#'
#' @examples
#' # Example usage with a contingency table:
#' pair_table <- table(c(1, 2, 2, 3), c(1, 1, 2, 2))
#' gen_number_max(pair_table)
#'
#' @export
gen_number_max <- function(x) {
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

  # Step 4: Calculate deltas for adjusting the values of selected cells
  delta1 <- min(table[RR[1], CC[1]], table[RR[2], CC[2]])
  delta2 <- -min(table[RR[1], CC[2]], table[RR[2], CC[1]])

  # Step 5: Create two modified tables with adjusted values
  table1 <- table
  table2 <- table

  # Modify table1
  table1[RR[1], CC[1]] <- table1[RR[1], CC[1]] - delta1
  table1[RR[2], CC[2]] <- table1[RR[2], CC[2]] - delta1
  table1[RR[1], CC[2]] <- table1[RR[1], CC[2]] + delta1
  table1[RR[2], CC[1]] <- table1[RR[2], CC[1]] + delta1

  # Modify table2
  table2[RR[1], CC[1]] <- table2[RR[1], CC[1]] - delta2
  table2[RR[2], CC[2]] <- table2[RR[2], CC[2]] - delta2
  table2[RR[1], CC[2]] <- table2[RR[1], CC[2]] + delta2
  table2[RR[2], CC[1]] <- table2[RR[2], CC[1]] + delta2

  # Step 6: Calculate mutual information for both modified tables
  mutinf1 <- get_mutual_information(table1)
  mutinf2 <- get_mutual_information(table2)

  # Step 7: Return the table with higher mutual information
  if (mutinf1 > mutinf2) {
    return(table1)
  } else {
    return(table2)
  }
}
