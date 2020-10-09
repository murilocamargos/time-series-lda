rm(list = ls())
source("src/tslda.r")

# Create a time series with an abrupt changepoint
y <- c(array(1, 100), array(5, 100))
y <- y + rnorm(200, 0, 1)

# Run TSLDA with an approximated number of clusters (3)
g <- tslda(y, 3)