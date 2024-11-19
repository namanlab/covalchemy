#' @title Objective Function for Structural Learning (SL)
#'
#' @description This function calculates the objective function for a structural learning task. It computes multiple components
#' such as the total variation (TV) distance between original and generated datasets (`X` vs. `X_prime`, `Y` vs. `Y_prime`),
#' the changes in regression coefficients (`beta0` and `beta1`), the R² differences for each category in `Z`, and the
#' inter-cluster centroid distances. The loss function combines these components using penalty parameters (`lambda1`, `lambda2`,
#' `lambda3`, `lambda4`).
#'
#' @param X A numeric vector representing the original values of `X`.
#' @param Y A numeric vector representing the original values of `Y`.
#' @param Z A categorical vector representing the categories for each observation.
#' @param X_prime A numeric vector representing the generated values of `X`.
#' @param Y_prime A numeric vector representing the generated values of `Y`.
#' @param p A numeric vector representing the target correlation values for each category in `Z`.
#' @param beta0_orig The original intercept value for the regression model.
#' @param beta1_orig The original slope value for the regression model.
#' @param lambda1 Penalty parameters to control the importance of different loss components.
#' @param lambda2 Penalty parameters to control the importance of different loss components.
#' @param lambda3 Penalty parameters to control the importance of different loss components.
#' @param lambda4 Penalty parameters to control the importance of different loss components.
#' @param R2_orig The original R² value for the model (not used directly in the calculation but might be for reference).
#' @param printc A boolean flag to control printing of intermediate values for debugging.
#' @return A numeric value representing the total loss calculated by the objective function.
#'
#' @examples
#' # Test Case 1: Simple random data with normal distribution
#' set.seed(123)
#' X <- rnorm(100)
#' Y <- rnorm(100)
#' Z <- sample(1:3, 100, replace = TRUE)
#' X_prime <- rnorm(100)
#' Y_prime <- rnorm(100)
#' p <- c(0.5, 0.7, 0.9)
#' beta0_orig <- 0
#' beta1_orig <- 1
#' lambda1 <- lambda2 <- lambda3 <- lambda4 <- 1
#' R2_orig <- 0.9
#' loss <- objective_function_SL(X, Y, Z, X_prime, Y_prime, p, beta0_orig, beta1_orig,
#'                                lambda1, lambda2, lambda3, lambda4, R2_orig)
#' print(loss)
#'
#' # Test Case 2: Skewed data with different categories and a larger lambda for penalty
#' X <- rexp(100)
#' Y <- rpois(100, lambda = 2)
#' Z <- sample(1:4, 100, replace = TRUE)
#' X_prime <- rnorm(100)
#' Y_prime <- rnorm(100)
#' p <- c(0.3, 0.5, 0.8, 0.6)
#' beta0_orig <- 0.5
#' beta1_orig <- 1.5
#' lambda1 <- lambda2 <- lambda3 <- 0.5
#' lambda4 <- 2
#' R2_orig <- 0.85
#' loss <- objective_function_SL(X, Y, Z, X_prime, Y_prime, p, beta0_orig, beta1_orig,
#'                                lambda1, lambda2, lambda3, lambda4, R2_orig)
#' print(loss)
#'
#' @export
objective_function_SL <- function(X_prime, Y_prime, X, Y, Z, p, beta0_orig, beta1_orig,
                                  lambda1, lambda2, lambda3, lambda4, R2_orig, printc = F) {
  # Step 1: Compute total variation (TV) distance between original and generated data for X and Y
  TV_X <- calculate_tv_distance_empirical(X, X_prime)  # TV distance for X
  TV_Y <- calculate_tv_distance_empirical(Y, Y_prime)  # TV distance for Y

  # Step 2: Fit linear regression model for the entire dataset to compute new regression coefficients
  fit <- lm(Y_prime ~ X_prime)  # Fit regression model
  beta0_new <- coef(fit)[1]  # New intercept
  beta1_new <- coef(fit)[2]  # New slope

  # Step 3: Compute the squared changes in regression coefficients relative to the original values
  delta_beta0 <- ((beta0_new - beta0_orig) / beta0_orig) ^ 2
  delta_beta1 <- ((beta1_new - beta1_orig) / beta1_orig) ^ 2

  # Step 4: Compute R² differences for each category in Z
  R2_diff <- 0
  categories <- unique(Z)
  for (cat in categories) {
    X_cat <- X[Z == cat]
    Y_cat <- Y[Z == cat]
    X_prime_cat <- X_prime[Z == cat]
    Y_prime_cat <- Y_prime[Z == cat]

    # Calculate correlation (R²) for current category
    rho_curr <- cor(Y_prime_cat, X_prime_cat)
    R2_diff <- R2_diff + (rho_curr - p[cat]) ^ 2  # Squared difference from target correlation

    if (printc) {
      cat("rho inter:", cat, "act:", rho_curr, "expected:", p[cat], "sqdiff:", (rho_curr - p[cat])^2, "\n")
    }
  }
  R2_diff <- R2_diff / length(categories)  # Average R² difference across categories

  # Step 5: Compute inter-cluster (centroid) distances
  cluster_centroids <- lapply(categories, function(cat) {
    X_cat_prime <- X_prime[Z == cat]
    Y_cat_prime <- Y_prime[Z == cat]
    return(c(mean(X_cat_prime), mean(Y_cat_prime)))  # Compute centroid for each category
  })

  inter_cluster_dist <- 0
  num_pairs <- 0
  for (i in 1:(length(categories) - 1)) {
    for (j in (i + 1):length(categories)) {
      centroid_i <- cluster_centroids[[i]]
      centroid_j <- cluster_centroids[[j]]

      # Calculate Euclidean distance between centroids
      dist_ij <- sqrt(sum((centroid_i - centroid_j) ^ 2))
      inter_cluster_dist <- inter_cluster_dist + dist_ij
      num_pairs <- num_pairs + 1
    }
  }

  # Step 6: Average inverse of inter-cluster distances
  inter_cluster_dist <- inter_cluster_dist / num_pairs

  # Step 7: Normalize inter-cluster distance using quantiles
  quantile_low <- 0.05  # 5th percentile for min
  quantile_high <- 0.95 # 95th percentile for max

  range_X <- quantile(X_prime, probs = c(quantile_low, quantile_high))
  range_Y <- quantile(Y_prime, probs = c(quantile_low, quantile_high))

  min_dist <- 0  # Minimum possible distance (perfect overlap)
  max_dist <- sqrt((range_X[2] - range_X[1])^2 + (range_Y[2] - range_Y[1])^2)  # Max possible distance

  # Step 8: Perform quantile-based min-max normalization
  normalized_inter_cluster_dist <- (inter_cluster_dist - min_dist) / (max_dist - min_dist)

  # Ensure that the normalized inter-cluster distance is between 0 and 1
  inter_cluster_dist <- max(0, min(normalized_inter_cluster_dist, 1))

  # Step 9: Calculate final loss function combining all components
  loss <- TV_X + TV_Y + lambda1 * delta_beta0^2 + lambda2 * delta_beta1^2 + lambda3 * R2_diff + lambda4 * inter_cluster_dist

  return(loss)  # Return the total loss value
}
