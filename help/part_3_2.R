log_odds_dc <- function(tab) {
  res <- ConDisPairs(tab)[c("C", "D")]
  concordant <- res$C + 1
  discordant <- res$D + 1
  if (discordant == 0) return(Inf)  # Handle division by zero
  return(log(concordant / discordant))
}


plot_log_odds <- function(matrices, new_matrices, names_matrices, log_odds_general = log_odds_dc) {
  before_log_odds <- sapply(matrices, log_odds_general)
  after_log_odds <- sapply(new_matrices, log_odds_general)
  names_matrices <- c(names_matrices, "Overall")
  overall_before <- log_odds_general(Reduce("+", matrices))
  overall_after <- log_odds_general(Reduce("+", new_matrices))
  before_log_odds <- c(before_log_odds, overall_before)
  after_log_odds <- c(after_log_odds, overall_after)
  data <- data.frame(matrix_name = c(names_matrices, names_matrices),
                     log_odds = c(before_log_odds, after_log_odds),
                     type = c(rep("Before", length(names_matrices)), rep("After", length(names_matrices))))
  print(data)
  plot_gen = data %>%
    mutate(type = factor(type, levels = c("Before", "After"), ordered = TRUE)) %>%
    ggplot(aes(x = matrix_name, y = log_odds, fill = type)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    xlab("Matrix") +
    ylab("Log Odds") +
    ggtitle("Log Odds Before and After Transformation")
  print(plot_gen)
}


softmax <- function(x) {
  # If any value is Inf, assign it a probability of 1, others 0
  if (any(is.infinite(x) & x > 0)) {
    probs <- rep(0, length(x))
    probs[which(x == Inf)] <- 1
    return(probs)
  }
  # Handle the case where all values are -Inf
  if (all(x == -Inf)) {
    return(rep(0, length(x)))
  }
  # Standard softmax calculation with numerical stability
  exp_x <- exp(x - max(x, na.rm = TRUE))  # Stability adjustment
  res = exp_x / sum(exp_x, na.rm = TRUE)
  return(res/sum(res))
}


augment_matrix_random_block <- function(table, delta) {
  nrows <- nrow(table)
  ncols <- ncol(table)
  iters <- 0
  repeat {
    iters <- iters + 1
    RR <- sort(sample(1:nrows, 2))
    CC <- sort(sample(1:ncols, 2))
    if (delta < 0 && table[RR[1], CC[1]] > 0 && table[RR[2], CC[2]] > 0) break
    if (delta > 0 && table[RR[1], CC[2]] > 0 && table[RR[2], CC[1]] > 0) break
    if (iters > 100) {
      delta <- 0
      break
    }
  }
  table[RR[1], CC[1]] <- table[RR[1], CC[1]] + delta
  table[RR[2], CC[2]] <- table[RR[2], CC[2]] + delta
  table[RR[1], CC[2]] <- table[RR[1], CC[2]] - delta
  table[RR[2], CC[1]] <- table[RR[2], CC[1]] - delta
  return(table)
}


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



set.seed(42)
matrices <- list(
  ta = matrix(c(512, 89, 313, 19), ncol = 2, byrow = TRUE),
  tb = matrix(c(353, 17, 207, 8), ncol = 2, byrow = TRUE),
  tc = matrix(c(120, 202, 205, 391), ncol = 2, byrow = TRUE),
  td = matrix(c(138, 131, 279, 244), ncol = 2, byrow = TRUE),
  te = matrix(c(53, 94, 138, 299), ncol = 2, byrow = TRUE),
  tf = matrix(c(22, 24, 351, 317), ncol = 2, byrow = TRUE)
)
levels_z <- names(matrices)
df_list <- lapply(seq_along(matrices), function(i) {
  mat <- matrices[[i]]
  z_level <- levels_z[i]
  df <- as.data.frame(as.table(mat))
  colnames(df) <- c("x", "y", "Freq")
  df$z <- z_level
  return(df)
})
# Combine all the data frames into one
final_df <- do.call(rbind, df_list)
expanded_df <- final_df[rep(1:nrow(final_df), final_df$Freq), c("x", "y", "z")]
result <- get_simpsons_paradox_d(expanded_df$x, expanded_df$y, expanded_df$z,
                                 manual_vec = c(-1, -1, -1, -1, -1, -1),
                                 target_overall = +1,
                                 margin = 0.2,  margin_overall = 0.2, max_n = 200)
table(expanded_df$x) - table(result$final_df$x)
table(expanded_df$y) - table(result$final_df$y)
table(expanded_df$z) - table(result$final_df$z)
