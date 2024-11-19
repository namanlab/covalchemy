#' Simpson's Paradox Transformation with Copula and Simulated Annealing
#'
#' This function simulates the Simpson's Paradox phenomenon by transforming data using Gaussian copulas,
#' optimizing the transformation with simulated annealing, and comparing the results.
#'
#' @param x A numeric vector of data points for variable X.
#' @param y A numeric vector of data points for variable Y.
#' @param z A categorical variable representing groups (e.g., factor or character vector).
#' @param corr_vector A vector of correlations for each category of z.
#' @param inv_cdf_type Type of inverse CDF transformation ("quantile_1", "quantile_4", "quantile_7", "quantile_8", "linear", "akima", "poly"). Default is "quantile_7".
#' @param sd_x Standard deviation for perturbations on X (default is 0.05).
#' @param sd_y Standard deviation for perturbations on Y (default is 0.05).
#' @param lambda1 Regularization parameter for simulated annealing (default is 1).
#' @param lambda2 Regularization parameter for simulated annealing (default is 1).
#' @param lambda3 Regularization parameter for simulated annealing (default is 1).
#' @param lambda4 Regularization parameter for simulated annealing (default is 1).
#' @param max_iter Maximum iterations for simulated annealing (default is 1000).
#' @param initial_temp Initial temperature for simulated annealing (default is 1.0).
#' @param cooling_rate Cooling rate for simulated annealing (default is 0.99).
#' @param order_vec Manual ordering of grids (default is NA, calculated automatically if not specified).
#'
#' @return A list containing:
#'   \item{df_all}{The final dataset with original, transformed, and annealed data.}
#'   \item{df_res}{A simplified version with only the optimized data.}
#'
#' @examples
#' set.seed(123)
#' n <- 300
#' z <- sample(c("A", "B", "C"), prob = c(0.3, 0.4, 0.3), size = n, replace = TRUE)
#' x <- rnorm(n, 10, sd = 5) + 5 * rbeta(n, 5, 3)
#' y <- 2 * x + rnorm(n, 5, sd = 4)
#' t <- c(-0.8, 0.8, -0.8)
#' res <- get_simpsons_paradox_c(x, y, z, t, sd_x = 0.07, sd_y = 0.07, lambda4 = 5)
#'
#' @export
get_simpsons_paradox_c <- function(x, y, z,
                                   corr_vector,   # Vector of correlations
                                   inv_cdf_type = "quantile_7",     # Inverse CDF function to transform samples
                                   sd_x = 0.05, sd_y = 0.05,   # Standard deviations for perturbations
                                   lambda1 = 1, lambda2 = 1, lambda3 = 1, lambda4 = 1,
                                   max_iter = 1000,            # Maximum iterations for simulated annealing
                                   initial_temp = 1.0,         # Initial temperature for annealing
                                   cooling_rate = 0.99,        # Cooling rate for annealing
                                   order_vec = NA) {           # to manually specify ordering of grids

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

  if (any(is.na(order_vec))) {
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
    labs(title = "Simpson's Paradox Copula Transformation",
         caption = str_c("Overall Correlation: ", round(cor(df_composition$x, df_composition$y), 3)))
  p2 <- ggMarginal(p_copula, type = "density")

  p_anneal <- ggplot(df_annealing) +
    geom_point(aes(x = x_optimized, y = y_optimized, color = z)) +
    theme_bw() + guides(color = "none") +
    labs(title = "Simulated Annealing Optimization",
         caption = str_c("Overall Correlation: ", round(cor(df_annealing$x_optimized, df_annealing$y_optimized), 3)))
  p3 <- ggMarginal(p_anneal, type = "density")

  # Plots:
  grid.arrange(p1, p2, p3, nrow = 1)
  res = df_final %>% select(x_optimized, y_optimized)

  # Return the combined dataframe containing all the transformations
  return(list(df_all = df_final, df_res = res))
}
