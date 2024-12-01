% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/objective_function_SL.R
\name{objective_function_SL}
\alias{objective_function_SL}
\title{Objective Function for Structural Learning (SL)}
\usage{
objective_function_SL(
  X_prime,
  Y_prime,
  X,
  Y,
  Z,
  p,
  beta0_orig,
  beta1_orig,
  lambda1,
  lambda2,
  lambda3,
  lambda4,
  R2_orig,
  printc = FALSE
)
}
\arguments{
\item{X_prime}{A numeric vector representing the generated values of \code{X}.}

\item{Y_prime}{A numeric vector representing the generated values of \code{Y}.}

\item{X}{A numeric vector representing the original values of \code{X}.}

\item{Y}{A numeric vector representing the original values of \code{Y}.}

\item{Z}{A categorical vector representing the categories for each observation.}

\item{p}{A numeric vector representing the target correlation values for each category in \code{Z}.}

\item{beta0_orig}{The original intercept value for the regression model.}

\item{beta1_orig}{The original slope value for the regression model.}

\item{lambda1}{Penalty parameters to control the importance of different loss components.}

\item{lambda2}{Penalty parameters to control the importance of different loss components.}

\item{lambda3}{Penalty parameters to control the importance of different loss components.}

\item{lambda4}{Penalty parameters to control the importance of different loss components.}

\item{R2_orig}{The original R² value for the model (not used directly in the calculation but might be for reference).}

\item{printc}{A boolean flag to control printing of intermediate values for debugging.}
}
\value{
A numeric value representing the total loss calculated by the objective function.
}
\description{
This function calculates the objective function for a structural learning task. It computes multiple components
such as the total variation (TV) distance between original and generated datasets (\code{X} vs. \code{X_prime}, \code{Y} vs. \code{Y_prime}),
the changes in regression coefficients (\code{beta0} and \code{beta1}), the R² differences for each category in \code{Z}, and the
inter-cluster centroid distances. The loss function combines these components using penalty parameters (\code{lambda1}, \code{lambda2},
\code{lambda3}, \code{lambda4}).
}
\examples{
# Test Case 1: Simple random data with normal distribution
set.seed(123)
X <- rnorm(100)
Y <- rnorm(100)
Z <- sample(1:3, 100, replace = TRUE)
X_prime <- rnorm(100)
Y_prime <- rnorm(100)
p <- c(0.5, 0.7, 0.9)
beta0_orig <- 0
beta1_orig <- 1
lambda1 <- lambda2 <- lambda3 <- lambda4 <- 1
R2_orig <- 0.9
loss <- objective_function_SL(X, Y, Z, X_prime, Y_prime, p, beta0_orig, beta1_orig,
                               lambda1, lambda2, lambda3, lambda4, R2_orig)
print(loss)

# Test Case 2: Skewed data with different categories and a larger lambda for penalty
X <- rexp(100)
Y <- rpois(100, lambda = 2)
Z <- sample(1:4, 100, replace = TRUE)
X_prime <- rnorm(100)
Y_prime <- rnorm(100)
p <- c(0.3, 0.5, 0.8, 0.6)
beta0_orig <- 0.5
beta1_orig <- 1.5
lambda1 <- lambda2 <- lambda3 <- 0.5
lambda4 <- 2
R2_orig <- 0.85
loss <- objective_function_SL(X, Y, Z, X_prime, Y_prime, p, beta0_orig, beta1_orig,
                               lambda1, lambda2, lambda3, lambda4, R2_orig)
print(loss)

}
