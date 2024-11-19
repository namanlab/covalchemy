# Function to calculate entropy of a pair
entropy_pair <- function(table) {
  total <- sum(table)
  probs <- table / total
  return(sum(-probs * log2(probs), na.rm = TRUE))
}

# Function to calculate mutual information
get_mutual_information <- function(table) {
  rowsums <- rowSums(table)
  colsums <- colSums(table)
  ent_x <- sum(-rowsums/sum(rowsums) * log2(rowsums/sum(rowsums)), na.rm = TRUE)
  ent_y <- sum(-colsums/sum(colsums) * log2(colsums/sum(colsums)), na.rm = TRUE)
  return(ent_x + ent_y - entropy_pair(table))
}

# Function to generate a new number for maximizing mutual information
gen_number_max <- function(x) {
  table <- x
  nrows <- nrow(table)
  ncols <- ncol(table)
  repeat {
    RR <- sample(1:nrows, 2)
    CC <- sample(1:ncols, 2)
    if (table[RR[1], CC[1]] > 0 && table[RR[2], CC[2]] > 0) break
  }

  delta1 <- min(table[RR[1], CC[1]], table[RR[2], CC[2]])
  delta2 <- -min(table[RR[1], CC[2]], table[RR[2], CC[1]])

  table1 <- table
  table2 <- table

  table1[RR[1], CC[1]] <- table1[RR[1], CC[1]] - delta1
  table1[RR[2], CC[2]] <- table1[RR[2], CC[2]] - delta1
  table1[RR[1], CC[2]] <- table1[RR[1], CC[2]] + delta1
  table1[RR[2], CC[1]] <- table1[RR[2], CC[1]] + delta1

  table2[RR[1], CC[1]] <- table2[RR[1], CC[1]] - delta2
  table2[RR[2], CC[2]] <- table2[RR[2], CC[2]] - delta2
  table2[RR[1], CC[2]] <- table2[RR[1], CC[2]] + delta2
  table2[RR[2], CC[1]] <- table2[RR[2], CC[1]] + delta2

  mutinf1 <- get_mutual_information(table1)
  mutinf2 <- get_mutual_information(table2)

  if (mutinf1 > mutinf2) {
    return(table1)
  } else {
    return(table2)
  }
}

# Function to generate a new number for minimizing mutual information
gen_number_min <- function(x) {
  table <- x
  nrows <- nrow(table)
  ncols <- ncol(table)
  repeat {
    RR <- sample(1:nrows, 2)
    CC <- sample(1:ncols, 2)
    if (table[RR[1], CC[1]] > 0 && table[RR[2], CC[2]] > 0) break
  }

  S <- (table[RR[1], CC[1]] + table[RR[2], CC[2]] + table[RR[1], CC[2]] + table[RR[2], CC[1]])
  delta <- round((table[RR[1], CC[1]]*table[RR[2], CC[2]] - table[RR[1], CC[2]]*table[RR[2], CC[1]]) / S)
  table[RR[1], CC[1]] <- table[RR[1], CC[1]] - delta
  table[RR[2], CC[2]] <- table[RR[2], CC[2]] - delta
  table[RR[1], CC[2]] <- table[RR[1], CC[2]] + delta
  table[RR[2], CC[1]] <- table[RR[2], CC[1]] + delta

  return(table)
}

# Function to generate a new number for stepwise modification
gen_number_1 <- function(x) {
  table <- x
  nrows <- nrow(table)
  ncols <- ncol(table)
  repeat {
    RR <- sample(1:nrows, 2)
    CC <- sample(1:ncols, 2)
    if (table[RR[1], CC[1]] > 0 && table[RR[2], CC[2]] > 0) break
  }

  delta <- 1

  table[RR[1], CC[1]] <- table[RR[1], CC[1]] - delta
  table[RR[2], CC[2]] <- table[RR[2], CC[2]] - delta
  table[RR[1], CC[2]] <- table[RR[1], CC[2]] + delta
  table[RR[2], CC[1]] <- table[RR[2], CC[1]] + delta

  return(table)
}

# Function to perform simulated annealing with stopping condition for target entropy
simulated_annealing_MI <- function(initial_table, obj, gen_fn, target,
                         max_n = 5000, temp = 10, maxim = T, readj = F) {
  best <- initial_table
  best_eval <- obj(best)
  curr <- best
  curr_eval <- best_eval
  n <- 0
  mutual_info_history <- data.frame(iteration = 0, mutual_info = curr_eval,
                                    type = ifelse(readj, "Readjusting", ifelse(maxim, "Maximizing", "Minimizing")))

  while (n < max_n) {
    if ((maxim && curr_eval >= target) || (!maxim && curr_eval <= target)) {
      print(paste("Target entropy reached:", curr_eval))
      break
    }

    cand <- gen_fn(curr)
    cand_eval <- obj(cand)
    if (maxim){
      if (cand_eval > best_eval) {
        best <- cand
        best_eval <- cand_eval
      }
      diff <- cand_eval - curr_eval
    } else {
      if (cand_eval < best_eval) {
        best <- cand
        best_eval <- cand_eval
      }
      diff <- -cand_eval + curr_eval
    }

    temp <- 0.9*temp #/ (n + 1)
    metropolis <- exp(diff / temp)
    if (diff > 0 || runif(1) < metropolis) {
      curr <- cand
      curr_eval <- cand_eval
    }
    n <- n + 1
    mutual_info_history <- rbind(mutual_info_history, data.frame(iteration = n, mutual_info = curr_eval,
                                                                 type =  ifelse(readj, "Readjusting", ifelse(maxim, "Maximizing", "Minimizing"))))
  }
  return(list(best, best_eval, n, mutual_info_history))
}

# Sinkhorn algorithm function
sinkhorn_algorithm <- function(initial_table, obj, max_iter = 500, tolerance = 1e-5) {
  S <- sum(initial_table)
  S_r <- rowSums(initial_table)
  S_c <- colSums(initial_table)
  n <- nrow(initial_table)
  m <- ncol(initial_table)

  # Initialize alpha and beta
  alpha <- rep(1, m)
  beta <- rep(1, n)

  curr_eval <- obj(initial_table)
  mutual_info_history <- data.frame(iteration = 0, mutual_info = curr_eval,
                                    type = "Minimizing")

  for (iter in 1:max_iter) {

    ones_mat <- matrix(1, nrow = m, ncol = n)

    # Update alpha
    alpha <- (S_c / (ones_mat %*% beta))
    # Update beta
    beta <- (S_r / (t(ones_mat) %*% alpha))

    # Convergence check
    updated_table <- t(diag(as.vector(alpha)) %*% ones_mat %*% diag(as.vector(beta)))
    curr_eval <- obj(updated_table)

    mutual_info_history <- rbind(mutual_info_history, data.frame(iteration = iter, mutual_info = curr_eval,
                                                                 type = "Minimizing"))
  }
  updated_table <- t(diag(as.vector(alpha)) %*% ones_mat %*% diag(as.vector(beta)))
  new_mut <- obj(updated_table)
  return(list(updated_table, new_mut, iter, mutual_info_history))
}


# Function to solve entropy by adjusting mutual information
get_target_entropy <- function(x, y, target_entropy, max_n = 10000, epsilon = 0.001) {
  table <- table(x, y)
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
  print(paste("Final Mutual Information:", final_mut))
  contingency_df <- as.data.frame(as.table(final_table))
  expanded_df <- contingency_df[rep(1:nrow(contingency_df), contingency_df$Freq), c("x", "y")]
  return(list(final_df = expanded_df, final_table = final_table, history = final_hist,
              max_mut = new_mut_max, min_mut = new_mut_min))
}

# Test cases:

# Increasing
set.seed(33)
df <- data.frame(
  x = sample(paste("Categ", 1:4), 10000, replace = TRUE),
  y = sample(paste("Categ", 10:4), 10000, replace = TRUE)
)

target_entropy <- 1 # Set your target entropy here
res <- get_target_entropy(df$x, df$y, target_entropy)
table(df$x)
table(res$final_df$x)
table(df$y)
table(res$final_df$y)
get_mutual_information(table(df$x, df$y))
get_mutual_information(table(res$final_df$x, res$final_df$y))


# Decreasing
set.seed(42)
library(stringr)
df <- data.frame(
  x = sample(str_c("Categ", 1:4), 10000, replace = TRUE)
) %>% mutate(y =
               ifelse(str_ends(x, "1"), sample(str_c("Categ", 10:8), 1, replace = TRUE),
                      ifelse(str_ends(x, "2"), sample(str_c("Categ", 8:6), 1, replace = TRUE),
                             sample(str_c("Categ", 5:4), 1, replace = TRUE))
               )
)
target_entropy <- 0.2
res <- get_target_entropy(df$x, df$y, target_entropy)
table(df$x)
table(res$final_df$x)
table(df$y)
table(res$final_df$y)
get_mutual_information(table(df$x, df$y))
get_mutual_information(table(res$final_df$x, res$final_df$y))

