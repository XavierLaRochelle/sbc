generate_y <- function(mu, kappa, j) {
  y <- c()
  for (n in 1:length(mu)) {
    p <- pnorm(mu[n] + kappa[j[n], ])
    theta <- c()
    theta[1] <- p[1]
    theta[2] <- p[2] - p[1]
    theta[3] <- 1 - p[2]
    y[n] <- sample(1:3, size = 1, prob = theta)
  }
  return(y)
}

source("./generate_params.r")

source("./generate_data.r")
