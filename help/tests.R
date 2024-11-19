################################################################################
################################## Part 1 ######################################
################################################################################

x1 = ChickWeight$weight
x2 = ChickWeight$Time
cor(x1, x2, method = "kendall")
res = get_target_corr(x1, x2, target_corr_kendall = 0,
                      copula_type = "t", inv_cdf_type = "poly", degree = 5)
cor(res$x1, res$x2, method = "kendall")
res = get_target_corr(x1, x2, target_corr_kendall = 0.99)
cor(res$x1, res$x2, method = "kendall")


################################################################################
################################## Part 2 ######################################
################################################################################

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




################################################################################
################################# Part 3.1 #####################################
################################################################################

# Generate data for two categories of z
# overall positive corr, works nice
set.seed(123)
n <- 300
z <- sample(c("A", "B", "C"), prob = c(0.3, 0.4, 0.3), size = n, replace = TRUE)
x <- rnorm(n, 10, sd = 5) + 5*rbeta(n, 5, 3)
y <- 2*x + rnorm(n, 5, sd = 4)
t = c(-0.8, 0.8, -0.8)
res = get_simpsons_paradox_c(x, y, z, t, sd_x = 0.07, sd_y = 0.07, lambda4 = 5)

################################################################################
################################# Part 3.2 #####################################
################################################################################


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

