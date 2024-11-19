#' Generate a New Number for Stepwise Modification
#'
#' This function modifies a given contingency table by swapping values between two cells
#' in a stepwise manner, where the change is fixed at a `delta` value of 1. The function
#' randomly selects two cells from the table and adjusts their values by subtracting and
#' adding the `delta` value.
#'
#' @param x A contingency table (numeric matrix or table).
#' @return A modified contingency table with stepwise adjustments.
#' @details
#' This function performs the following steps:
#' 1. Randomly selects two rows and two columns from the table.
#' 2. Ensures that the selected cells have non-zero values.
#' 3. Adjusts the values of the selected cells by subtracting 1 from two cells and
#'    adding 1 to the other two.
#' 4. Returns the modified table with stepwise adjustments.
#'
#' @examples
#' # Example usage with a contingency table:
#' pair_table <- table(c(1, 2, 2, 3), c(1, 1, 2, 2))
#' gen_number_1(pair_table)
#'
#' @export
gen_number_1 <- function(x) {
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

  # Step 4: Set delta value for modification
  delta <- 1

  # Step 5: Modify the values of the selected cells
  table[RR[1], CC[1]] <- table[RR[1], CC[1]] - delta
  table[RR[2], CC[2]] <- table[RR[2], CC[2]] - delta
  table[RR[1], CC[2]] <- table[RR[1], CC[2]] + delta
  table[RR[2], CC[1]] <- table[RR[2], CC[1]] + delta

  # Step 6: Return the modified table
  return(table)
}
