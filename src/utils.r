generate_y <- function(mu, kappa, j) {
  y <- c()
  for (n in seq_along(mu)) {
    p <- pnorm(mu[n] + kappa[j[n], ])
    theta <- c()
    theta[1] <- p[1]
    theta[2] <- p[2] - p[1]
    theta[3] <- 1 - p[2]
    y[n] <- sample(1:3, size = 1, prob = theta)
  }
  return(y)
}

mo <- function(x, scale, j) {
  K <- ncol(x)
  N <- nrow(x)
  eta <- matrix(NA, nrow = N, ncol = K)
  for (k in 1:K) {
    for (n in 1:N) {
      if (x[n, k] == 1) {
        eta[n, k] <- 0
      } else if (x[n, k] == 3) {
        eta[n, k] <- 1
      } else {
        eta[n, k] <- scale[j[n], k, 1]
      }
    }
  }
  return(eta)
}
