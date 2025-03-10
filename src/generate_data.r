generate_data <- function(conf, params) {
  set.seed(conf$seed)
  N <- conf$N
  J <- conf$J
  K <- conf$K

  d <- list()
  j <- sample(1:J, size = N, replace = TRUE)
  d$j <- j
  X <- W <- matrix(NA, nrow = N, ncol = K)
  for (k in 1:K) {
    X[, k] <- sample(1:3, size = N, replace = TRUE)
    W[, k] <- sample(1:3, size = N, replace = TRUE)
  }
  d$X <- X
  d$W <- W

  kappa <- params$kappa
  beta <- params$beta
  gamma <- params$gamma
  lambda <- params$lambda
  omega <- params$omega
  delta <- params$delta
  zeta <- params$zeta

  mu <- c()
  eta_X <- mo(X, delta, j)
  eta_W <- mo(W, zeta, j)
  for (n in 1:N) {
    mu[n] <- eta_X[n, ] %*% (beta + gamma[j[n], ]) +
      (eta_X[n, ] * eta_W[n, ]) %*% (lambda + omega[j[n], ])
  }
  d$y <- generate_y(mu, kappa, j)

  return(d)
}
