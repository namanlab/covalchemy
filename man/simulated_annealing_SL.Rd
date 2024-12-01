% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulated_annealing_SL.R
\name{simulated_annealing_SL}
\alias{simulated_annealing_SL}
\title{Simulated Annealing Optimization with Categorical Variable and R^2 Differences}
\usage{
simulated_annealing_SL(
  X,
  Y,
  Z,
  X_st,
  Y_st,
  p,
  sd_x = 0.05,
  sd_y = 0.05,
  lambda1 = 1,
  lambda2 = 1,
  lambda3 = 1,
  lambda4 = 1,
  max_iter = 1000,
  initial_temp = 1,
  cooling_rate = 0.99
)
}
\arguments{
\item{X}{A numeric vector or matrix of input data (independent variable).}

\item{Y}{A numeric vector of the dependent variable (target).}

\item{Z}{A categorical variable (vector), used for grouping data in the analysis.}

\item{X_st}{A numeric vector of starting values for the composition method of X.}

\item{Y_st}{A numeric vector of starting values for the composition method of Y.}

\item{p}{A numeric vector representing the target R^2 values for each category in Z.}

\item{sd_x}{Standard deviation for the noise added to X during the perturbation (default is 0.05).}

\item{sd_y}{Standard deviation for the noise added to Y during the perturbation (default is 0.05).}

\item{lambda1}{Regularization parameter for the total variation distance term (default is 1).}

\item{lambda2}{Regularization parameter for the coefficient difference term (default is 1).}

\item{lambda3}{Regularization parameter for the R^2 difference term (default is 1).}

\item{lambda4}{Regularization parameter for the inter-cluster distance term (default is 1).}

\item{max_iter}{Maximum number of iterations for the annealing process (default is 1000).}

\item{initial_temp}{Initial temperature for the annealing process (default is 1.0).}

\item{cooling_rate}{The rate at which the temperature cools down during annealing (default is 0.99).}
}
\value{
A list with the optimized values of X_prime and Y_prime.
}
\description{
This function implements the Simulated Annealing algorithm to optimize a solution
based on the total variation distance, changes in regression coefficients,
R-squared differences, and inter-cluster distance, with respect to a set of
categorical and continuous variables.
}
