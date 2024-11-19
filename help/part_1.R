generate_gaussian_copula_samples <- function(n, d, rho_matrix) {
  # Step 1: Generate multivariate normal samples
  mean_vector <- rep(0, d)  # Mean vector for multivariate normal
  mvn_samples <- rmvnorm(n = n, mean = mean_vector, sigma = rho_matrix)

  # Step 2: Transform samples to uniform using the CDF of the standard normal
  uniform_samples <- pnorm(mvn_samples)

  return(uniform_samples)
}

gaussian_copula_two_vars <- function(n, p) {
  rho_matrix <- matrix(c(1, p, p, 1), nrow = 2, ncol = 2)
  smpls <- generate_gaussian_copula_samples(n, 2, rho_matrix)
  return(smpls)
}

generate_t_copula_samples <- function(n, d, rho_matrix, df) {
  # Step 1: Generate multivariate t samples
  mvn_samples <- rmvt(n = n, sigma = rho_matrix, df = df)

  # Step 2: Transform samples to uniform using the CDF of the t-distribution
  uniform_samples <- pt(mvn_samples, df = df)

  return(uniform_samples)
}

t_copula_two_vars <- function(n, p) {
  rho_matrix <- matrix(c(1, p, p, 1), nrow = 2, ncol = 2)
  smpls <- generate_t_copula_samples(n, 2, rho_matrix, 5)
  return(smpls)
}


# Inverse CDF generator using quantiles
genCDFInv_quantile <- function(data, type = 1) {
  return(function(p) {quantile(data, p, type = type)})
}

# Inverse CDF generator using linear interpolation
genCDFInv_linear <- function(X) {
  ecdf1 <- ecdf(X)
  U <- sort((ecdf1(X)))
  Finv <- sort(X)
  ak2 <- approxfun(U, Finv, method = "linear", rule = 2)
  return(ak2)
}

# Inverse CDF generator using Akima spline interpolation
genCDFInv_akima <- function(X) {
  ecdf1 <- ecdf(X)
  U <- sort((ecdf1(X)))
  Finv <- sort(X)
  ak2 <- function(x) {aspline(U, Finv, x)$y}
  return(ak2)
}

# Inverse CDF generator using polynomial regression
genCDFInv_poly <- function(data, degree) {
  x <- sort((data))
  y <- ecdf(data)(x)
  poly_fit <- lm(x ~ poly(y, degree = degree, raw = TRUE))
  cdf_poly <- function(y) predict(poly_fit, newdata = data.frame(y = y))
  return(cdf_poly)
}


get_target_corr <- function(x1, x2, target_corr_kendall,
                            copula_type = "gaussian", inv_cdf_type = "quantile_7",
                            degree = 10) {
  n <- length(x1)

  # Generate copula samples
  if (copula_type == "gaussian") {
    target_rho <- sin(target_corr_kendall * pi / 2) # Greiner's equality
    res <- gaussian_copula_two_vars(n, target_rho)
  } else if (copula_type == "t") {
    target_rho <- sin(target_corr_kendall * pi / 2) # Greiner's equality
    res <- t_copula_two_vars(n, target_rho)
  } else {
    stop("Unsupported copula type")
  }

  # Apply chosen inverse CDF
  inv_cdf_d1 <- switch(inv_cdf_type,
                       "quantile_1" = genCDFInv_quantile(x1, 1), # stepwise
                       "quantile_4" = genCDFInv_quantile(x1, 4), # linear interpolation
                       "quantile_7" = genCDFInv_quantile(x1, 7), # default
                       "quantile_8" = genCDFInv_quantile(x1, 8), # median-unbiased
                       "linear" = genCDFInv_linear(x1),
                       "akima" = genCDFInv_akima(x1),
                       "poly" = genCDFInv_poly(x1, degree = degree))

  inv_cdf_d2 <- switch(inv_cdf_type,
                       "quantile_1" = genCDFInv_quantile(x2, 1), # stepwise
                       "quantile_4" = genCDFInv_quantile(x2, 4), # linear interpolation
                       "quantile_7" = genCDFInv_quantile(x2, 7), # default
                       "quantile_8" = genCDFInv_quantile(x2, 8), # median-unbiased
                       "linear" = genCDFInv_linear(x2),
                       "akima" = genCDFInv_akima(x2),
                       "poly" = genCDFInv_poly(x2, degree = degree))

  F1Inv <- Vectorize(inv_cdf_d1)
  F2Inv <- Vectorize(inv_cdf_d2)

  x1mod <- F1Inv(res[, 1])
  x2mod <- F2Inv(res[, 2])

  return(list(x1 = x1mod, x2 = x2mod))
}


# sample:
x1 = ChickWeight$weight
x2 = ChickWeight$Time
cor(x1, x2, method = "kendall")
res = get_target_corr(x1, x2, target_corr_kendall = 0,
                      copula_type = "t", inv_cdf_type = "poly", degree = 5)
cor(res$x1, res$x2, method = "kendall")
res = get_target_corr(x1, x2, target_corr_kendall = 0.99)
cor(res$x1, res$x2, method = "kendall")
