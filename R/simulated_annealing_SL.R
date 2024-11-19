#' Simulated Annealing Optimization with Categorical Variable and R^2 Differences
#'
#' This function implements the Simulated Annealing algorithm to optimize a solution
#' based on the total variation distance, changes in regression coefficients,
#' R-squared differences, and inter-cluster distance, with respect to a set of
#' categorical and continuous variables.
#'
#' @param X A numeric vector or matrix of input data (independent variable).
#' @param Y A numeric vector of the dependent variable (target).
#' @param Z A categorical variable (vector), used for grouping data in the analysis.
#' @param X_st A numeric vector of starting values for the composition method of X.
#' @param Y_st A numeric vector of starting values for the composition method of Y.
#' @param p A numeric vector representing the target R^2 values for each category in Z.
#' @param sd_x Standard deviation for the noise added to X during the perturbation (default is 0.05).
#' @param sd_y Standard deviation for the noise added to Y during the perturbation (default is 0.05).
#' @param lambda1 Regularization parameter for the total variation distance term (default is 1).
#' @param lambda2 Regularization parameter for the coefficient difference term (default is 1).
#' @param lambda3 Regularization parameter for the R^2 difference term (default is 1).
#' @param lambda4 Regularization parameter for the inter-cluster distance term (default is 1).
#' @param max_iter Maximum number of iterations for the annealing process (default is 1000).
#' @param initial_temp Initial temperature for the annealing process (default is 1.0).
#' @param cooling_rate The rate at which the temperature cools down during annealing (default is 0.99).
#'
#' @return A list with the optimized values of X_prime and Y_prime.
#' @export
simulated_annealing_SL <- function(X, Y, Z, X_st, Y_st, p, sd_x = 0.05, sd_y = 0.05,
                                lambda1 = 1, lambda2 = 1, lambda3 = 1, lambda4 = 1,
                                max_iter = 1000, initial_temp = 1.0, cooling_rate = 0.99) {

  # Initialize X_prime and Y_prime from the starting composition method output
  X_prime <- X_st
  Y_prime <- Y_st
  best_X_prime <- X_prime
  best_Y_prime <- Y_prime

  # Compute original regression coefficients (beta0 and beta1) using the full dataset
  fit_orig <- lm(Y ~ X)
  beta0_orig <- coef(fit_orig)[1]
  beta1_orig <- coef(fit_orig)[2]

  # Compute R-squared for each category in the categorical variable Z
  categories <- unique(Z)
  R2_orig <- sapply(categories, function(cat) {
    fit_cat <- lm(Y[Z == cat] ~ X[Z == cat])
    summary(fit_cat)$r.squared
  })

  # Initialize best_loss using the objective function with initial X_prime and Y_prime
  best_loss <- objective_function_SL(X_prime, Y_prime, X, Y, Z, p, beta0_orig, beta1_orig,
                                     lambda1, lambda2, lambda3, lambda4, R2_orig)

  # Start the simulated annealing process with the initial temperature
  temp <- initial_temp
  for (i in 1:max_iter) {
    # Generate new solutions by adding random noise to X_prime and Y_prime
    new_X_prime <- X_prime + rnorm(length(X_prime), 0, sd_x)
    new_Y_prime <- Y_prime + rnorm(length(Y_prime), 0, sd_y)

    # Compute the new loss using the updated values
    new_loss <- objective_function_SL(new_X_prime, new_Y_prime, X, Y, Z, p, beta0_orig,
                                      beta1_orig, lambda1, lambda2, lambda3, lambda4, R2_orig)

    # Accept the new solution based on the Metropolis criterion
    if (new_loss < best_loss || runif(1) < exp((best_loss - new_loss) / temp)) {
      X_prime <- new_X_prime
      Y_prime <- new_Y_prime
      best_loss <- new_loss
      best_X_prime <- X_prime
      best_Y_prime <- Y_prime
    }

    # Decrease the temperature according to the cooling rate
    temp <- temp * cooling_rate

    # Check if the temperature is low enough to stop
    if (temp < 1e-6) {
      break
    }
  }

  # Return the best solution found
  return(list(X_prime = best_X_prime, Y_prime = best_Y_prime))
}

