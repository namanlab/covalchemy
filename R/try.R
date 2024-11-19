install.packages(c("devtools", "roxygen2", "usethis", "testthat"))

devtools::document() # manual

devtools::check()


library(tidyverse)
library(mvtnorm)
library(akima)
library(clue)
library(ggExtra)
library(gridExtra)
library(DescTools)

devtools::document()
