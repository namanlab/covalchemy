# covalchemy

This package provides functions for manipulating data, including techniques for:

-   Modifying correlations for continuous variables: The `get_target_corr` function allows you to adjust the Kendall's tau correlation between two continuous variables using copula-based methods. You can specify the target correlation, copula type (Gaussian or t), and inverse CDF transformation method.
-   Controlling entropy in categorical variables: The `get_target_entropy` function helps you achieve a desired level of entropy (mutual information) between two categorical variables. It works by iteratively adjusting the contingency table using simulated annealing.
-   Simulating Simpson's Paradox with continuous variables: The `get_simpsons_paradox_c` function allows you to explore the Simpson's Paradox phenomenon in continuous data. It transforms data using Gaussian copulas and simulated annealing to create a scenario where the overall trend contradicts subgroup trends.
-   Simulating Simpson's Paradox-like effects in categorical data: The `get_simpsons_paradox_d` function enables you to modify contingency tables for categorical data to create or highlight a Simpson's Paradox-like effect. It employs simulated annealing to adjust log-odds values while respecting specific constraints.

## Installation

To install the package, you can use the devtools package:

``` r
# Install devtools package if not already available
if (!requireNamespace("devtools")) install.packages("devtools")
devtools::install_github("namanlab/covalchemy")  
```

Once installed, you can load the package and use its functions in your R scripts:

``` r
library(covalchemy)

# Example 1: Modifying correlation
x <- rnorm(100)
y <- rnorm(100)
target_corr <- 0.5
res <- get_target_corr(x, y, target_corr)
modified_x <- res$x1
modified_y <- res$x2

# Example 2: Controlling entropy
df <- data.frame(x = sample(c("A", "B", "C"), 1000, replace = TRUE),
                  y = sample(c("D", "E", "F"), 1000, replace = TRUE))
target_entropy <- 1.5
result <- get_target_entropy(df$x, df$y, target_entropy)
final_df <- result$final_df
```

## Additional notes

This readme provides a general overview of the package's functionalities. Refer to the function documentation within the package for detailed information on arguments, return values, and specific usage examples.

## Acknowledgments

This package was developed as part of the DSA42288S Final Year Project. I would like to express my gratitude to my supervisor, Dr. Vikneswaran Gopal, for his invaluable guidance, support, and mentorship throughout this project. I am also grateful to the faculty and staff at NUS for their continuous support, and to my family for their encouragement along the way.

## Contact

For any questions or inquiries, please contact me at `naman.agr03@gmail.com`.
