calculate_tv_distance_empirical <- function(original_data, generated_data) {
  # Create empirical CDFs for original and generated data
  original_ecdf <- ecdf(original_data)
  generated_ecdf <- ecdf(generated_data)
  # Define the grid over which to calculate the distance
  x_values <- sort(unique(c(original_data, generated_data)))
  # Calculate the TV distance
  tv_distance <- max(abs(original_ecdf(x_values) - generated_ecdf(x_values)))
  return(tv_distance)
}


get_optimal_grid <- function(x, y, z){
  # Calculate quantiles for x and y based on z distribution
  z_levels <- sort(unique(z))
  x_quantiles <- quantile(x, probs = cumsum(table(z) / length(z)))
  y_quantiles <- quantile(y, probs = cumsum(table(z) / length(z)))
  k = length(unique(z))
  res = matrix(rep(0, k^2), nrow = k)
  for (i in seq_along(z_levels)) {
    x_range <- if (i == 1) {
      which(x <= x_quantiles[i])
    } else {
      which(x > x_quantiles[i - 1] & x <= x_quantiles[i])
    }
    for (j in seq_along(z_levels)) {
      temp_y = y[x_range]
      y_range <- if (j == 1) {
        which(temp_y <= y_quantiles[j])
      } else {
        which(temp_y > y_quantiles[j - 1] & temp_y <= y_quantiles[j])
      }
      res[i, j] = length(y_range)

    }
  }
  # Convert the matrix to a cost matrix
  M <- max(res)
  cost_matrix <- M - res
  print(res)

  # Solve the assignment problem using the Hungarian algorithm
  assignment <- solve_LSAP(cost_matrix)

  # Return the row and column indices of the selected cells
  return(as.vector(assignment))
}


objective_function_SL <- function(X_prime, Y_prime, X, Y, Z, p, beta0_orig, beta1_orig,
                               lambda1, lambda2, lambda3, lambda4, R2_orig, printc = F) {
  # Compute total variation distance
  TV_X <- calculate_tv_distance_empirical(X, X_prime)
  TV_Y <- calculate_tv_distance_empirical(Y, Y_prime)

  # Fit regression model for the entire dataset
  fit <- lm(Y_prime ~ X_prime)
  beta0_new <- coef(fit)[1]
  beta1_new <- coef(fit)[2]

  # Compute the changes in regression coefficients
  delta_beta0 <- ((beta0_new - beta0_orig)/beta0_orig)^2
  delta_beta1 <- ((beta1_new - beta1_orig)/beta1_orig)^2

  # Compute R^2 for each category in Z
  R2_diff <- 0
  categories <- unique(Z)
  for (cat in categories) {
    X_cat <- X[Z == cat]
    Y_cat <- Y[Z == cat]
    X_prime_cat <- X_prime[Z == cat]
    Y_prime_cat <- Y_prime[Z == cat]

    # New R^2 for this category
    rho_curr <- cor(Y_prime_cat, X_prime_cat)
    R2_diff <- R2_diff + (rho_curr - p[cat])^2 # R2_orig_cat is target corr at start
    if (printc) {
      cat("rho inter:", cat, "act:", rho_curr, "expected:", p[cat], "sqdiff:", (rho_curr - p[cat])^2, "\n")
    }
  }
  R2_diff <- R2_diff / length(categories)

  # Compute inter-cluster (centroid) distances
  cluster_centroids <- lapply(categories, function(cat) {
    X_cat_prime <- X_prime[Z == cat]
    Y_cat_prime <- Y_prime[Z == cat]
    return(c(mean(X_cat_prime), mean(Y_cat_prime))) # Centroid for (X, Y) in category 'cat'
  })

  inter_cluster_dist <- 0
  num_pairs <- 0
  for (i in 1:(length(categories) - 1)) {
    for (j in (i + 1):length(categories)) {
      centroid_i <- cluster_centroids[[i]]
      centroid_j <- cluster_centroids[[j]]

      # Compute Euclidean distance between centroids
      dist_ij <- sqrt(sum((centroid_i - centroid_j)^2))
      inter_cluster_dist <- inter_cluster_dist + dist_ij
      num_pairs <- num_pairs + 1
    }
  }

  # Average inverse of inter-cluster distances
  inter_cluster_dist <- inter_cluster_dist / num_pairs

  # Use quantiles to compute a robust min and max distance (to handle outliers)
  quantile_low <- 0.05  # 5th percentile for minimum
  quantile_high <- 0.95 # 95th percentile for maximum

  # Quantile-based min and max distances for normalization
  range_X <- quantile(X_prime, probs = c(quantile_low, quantile_high))
  range_Y <- quantile(Y_prime, probs = c(quantile_low, quantile_high))

  # Min and max based on quantile ranges
  min_dist <- 0 # Since min_dist is still 0 (clusters perfectly overlap)
  max_dist <- sqrt((range_X[2] - range_X[1])^2 + (range_Y[2] - range_Y[1])^2) # Max possible distance based on quantiles

  # Perform quantile-based min-max normalization
  normalized_inter_cluster_dist <- (inter_cluster_dist - min_dist) / (max_dist - min_dist)

  # Ensure the normalized value is between 0 and 1
  inter_cluster_dist <- max(0, min(normalized_inter_cluster_dist, 1))

  # Loss function: combine TV distance, coefficient changes, R^2 differences, and cluster separation penalty
  loss <- TV_X + TV_Y + lambda1 * delta_beta0^2 + lambda2 * delta_beta1^2 + lambda3 * R2_diff + lambda4 * inter_cluster_dist
  return(loss)
}

# Simulated Annealing including categorical variable Z and R^2 differences
simulated_annealing_SL <- function(X, Y, Z, X_st, Y_st, p, sd_x = 0.05, sd_y = 0.05,
                                lambda1 = 1, lambda2 = 1, lambda3 = 1, lambda4 = 1,
                                max_iter = 1000, initial_temp = 1.0, cooling_rate = 0.99) {
  X_prime <- X_st # start form the composition method output
  Y_prime <- Y_st # start form the composition method output
  best_X_prime <- X_prime
  best_Y_prime <- Y_prime

  # Original regression coefficients
  fit_orig <- lm(Y ~ X)
  beta0_orig <- coef(fit_orig)[1]
  beta1_orig <- coef(fit_orig)[2]

  # Original R^2 for each category of Z
  categories <- unique(Z)
  R2_orig <- sapply(categories, function(cat) {
    fit_cat <- lm(Y[Z == cat] ~ X[Z == cat])
    summary(fit_cat)$r.squared
  })

  best_loss <- objective_function_SL(X_prime, Y_prime, X, Y, Z, p, beta0_orig, beta1_orig, lambda1, lambda2, lambda3, lambda4, R2_orig)

  temp <- initial_temp
  for (i in 1:max_iter) {
    # Generate a new solution (random perturbation)
    new_X_prime <- X_prime + rnorm(length(X_prime), 0, sd_x)
    new_Y_prime <- Y_prime + rnorm(length(Y_prime), 0, sd_y)

    # Compute the new loss
    new_loss <- objective_function_SL(new_X_prime, new_Y_prime, X, Y, Z, p, beta0_orig, beta1_orig, lambda1, lambda2, lambda3, lambda4, R2_orig)

    # Accept new solution based on probability
    if (new_loss < best_loss || runif(1) < exp((best_loss - new_loss) / temp)) {
      X_prime <- new_X_prime
      Y_prime <- new_Y_prime
      best_loss <- new_loss
      best_X_prime <- X_prime
      best_Y_prime <- Y_prime
    }

    # Cool down
    temp <- temp * cooling_rate

    # Convergence check
    if (temp < 1e-6) {
      break
    }
  }
  return(list(X_prime = best_X_prime, Y_prime = best_Y_prime))
}

################################################################################
############################### Combined Workflow ##############################
################################################################################


get_simpsons_paradox_c <- function(x, y, z,
                        corr_vector,   # Vector of correlations
                        inv_cdf_type = "quantile_7",     # Inverse CDF function to transform samples
                        sd_x = 0.05, sd_y = 0.05,   # Standard deviations for perturbations
                        lambda1 = 1, lambda2 = 1, lambda3 = 1, lambda4 = 1,
                        max_iter = 1000,            # Maximum iterations for simulated annealing
                        initial_temp = 1.0,         # Initial temperature for annealing
                        cooling_rate = 0.99,  # Cooling rate for annealing
                        order_vec = NA) {     # to mnaulayy psecify ordering of grids

  cat("Starting transformation process...\n")

  # Apply chosen inverse CDF
  FInvFunc <- switch(inv_cdf_type,
                       "quantile_1" = function(x1) {genCDFInv_quantile(x1, 1)}, # stepwise
                       "quantile_4" = function(x1) {genCDFInv_quantile(x1, 4)}, # linear interpolation
                       "quantile_7" = function(x1) {genCDFInv_quantile(x1, 7)}, # default
                       "quantile_8" = function(x1) {genCDFInv_quantile(x1, 8)}, # median-unbiased
                       "linear" = function(x1) {genCDFInv_linear(x1)},
                       "akima" = function(x1) {genCDFInv_akima(x1)},
                       "poly" = function(x1) {genCDFInv_poly(x1, degree = degree)})

  # Transform correlation vector for Gaussian copula
  p_corr_vec <- sin(corr_vector * pi / 2)

  # Get unique levels of z and sort data by z
  z_levels <- sort(unique(z))
  names(p_corr_vec) <- z_levels
  z <- sort(z)

  if (any(is.na(order_vec))){
    order_vec <- get_optimal_grid(x, y, z)
  }

  # Calculate quantiles for x and y based on z distribution
  x_quantiles <- quantile(x, probs = cumsum(table(z) / length(z)))
  y_quantiles <- quantile(y, probs = cumsum(table(z) / length(z)))

  df_composition <- data.frame()  # Initialize result dataframe

  # Iterate over each category of z
  for (i in seq_along(z_levels)) {
    z_cat <- z_levels[i]
    prob_val <- table(z)[i] / length(z)   # Proportion of category in z

    cat("Processing category:", z_cat, "with proportion:", prob_val, "\n")

    # Generate samples using the Gaussian copula
    p_corr <- p_corr_vec[i]
    n_samples <- sum(z == z_cat)  # Number of samples for current category
    copula_samples <- gaussian_copula_two_vars(n_samples, p_corr)  # Generate copula samples

    # Filter x and y samples for the current quantile range
    x_range <- if (i == 1) {
      which(x <= x_quantiles[i])
    } else {
      which(x > x_quantiles[i - 1] & x <= x_quantiles[i])
    }
    j = order_vec[i]
    y_range <- if (j == 1) {
      which(y <= y_quantiles[j])
    } else {
      which(y > y_quantiles[j - 1] & y <= y_quantiles[j])
    }

    x_samples <- x[x_range]
    y_samples <- y[y_range]

    # Generate inverse CDF functions for x and y samples
    F1Inv <- Vectorize(FInvFunc(x_samples))
    F2Inv <- Vectorize(FInvFunc(y_samples))

    # Transform copula samples using the inverse CDF functions
    x_transformed <- F1Inv(copula_samples[, 1])
    y_transformed <- F2Inv(copula_samples[, 2])

    # Append transformed data to the result dataframe
    df_composition <- rbind(df_composition, data.frame(x = x_transformed, y = y_transformed, z = z_cat))
  }

  cat("Piecewise copula transformation complete.\n")

  # Apply simulated annealing to further optimize the samples
  cat("Starting simulated annealing optimization...\n")
  res_anneal <- simulated_annealing_SL(x, y, z,
                                    df_composition$x, df_composition$y, p_corr_vec,
                                    sd_x, sd_y,
                                    lambda1, lambda2, lambda3, lambda4,
                                    max_iter, initial_temp, cooling_rate)

  # Create a final dataframe with the optimized x and y values
  df_annealing <- data.frame(x_optimized = res_anneal$X_prime,
                             y_optimized = res_anneal$Y_prime,
                             z = z)

  cat("Simulated annealing complete.\n")

  # Combine original, copula-transformed, and annealed data into a single dataframe
  df_final <- data.frame(x_original = x, y_original = y, z = z,
                         x_transformed = df_composition$x,
                         y_transformed = df_composition$y,
                         x_optimized = df_annealing$x_optimized,
                         y_optimized = df_annealing$y_optimized)

  # Calculate total variation (TV) distance before and after annealing
  tv_initial <- calculate_tv_distance_empirical(x, df_composition$x) +
    calculate_tv_distance_empirical(y, df_composition$y)
  tv_final <- calculate_tv_distance_empirical(x, df_annealing$x_optimized) +
    calculate_tv_distance_empirical(y, df_annealing$y_optimized)

  # Print TV distances
  cat("Total Variation Before Annealing: ", tv_initial, "\n")
  cat("Total Variation After Annealing: ", tv_final, "\n")


  # Plot original, copula-transformed, and simulated annealing results
  p_orig <- ggplot(data.frame(x = x, y = y, z = z)) +
    geom_point(aes(x = x, y = y, color = z)) +
    theme_bw() + guides(color = "none") +
    labs(title = "Original Data",
         caption = str_c("Overall Correlation: ", round(cor(x, y), 3)))
  p1 <- ggMarginal(p_orig, type = "density")

  p_copula <- ggplot(df_composition) +
    geom_point(aes(x = x, y = y, color = z)) +
    theme_bw() + guides(color = "none") +
    labs(title = "After Piecewise Copulas",
         caption = str_c("Overall Correlation: ",
                         round(cor(df_composition$x, df_composition$y), 3),
                         "\nTV Distance: ", round(tv_initial, 3)))
  p2 <- ggMarginal(p_copula, type = "density")

  p_anneal <- ggplot(df_annealing) +
    geom_point(aes(x = x_optimized, y = y_optimized, color = z)) +
    theme_bw() + guides(color = "none") +
    labs(x = "x", y = "y", title = "After Simulated Annealing",
         caption = str_c("Overall Correlation: ",
                         round(cor(df_annealing$x, df_annealing$y), 3),
                         "\nTV Distance: ", round(tv_final, 3)))
  p3 <- ggMarginal(p_anneal, type = "density")

  # Plots:
  grid.arrange(p1, p2, p3, nrow = 1)

  res = df_final %>% select(x_optimized, y_optimized)

  # Return the combined dataframe containing all the transformations
  return(list(df_all = df_final, df_res = res))
}



# Generate data for two categories of z
# overall positive corr, works nice
set.seed(123)
n <- 300
z <- sample(c("A", "B", "C"), prob = c(0.3, 0.4, 0.3), size = n, replace = TRUE)
x <- rnorm(n, 10, sd = 5) + 5*rbeta(n, 5, 3)
y <- 2*x + rnorm(n, 5, sd = 4)
t = c(-0.8, 0.8, -0.8)
res = get_simpsons_paradox_c(x, y, z, t, sd_x = 0.07, sd_y = 0.07, lambda4 = 5)
