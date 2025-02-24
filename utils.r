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

generate_params <- function(conf) {
  require(extraDistr)
  set.seed(conf$seed)
  J <- conf$J
  K <- conf$K
  mu_r2 <- conf$prior$mu_r2
  phi_r2 <- conf$prior$phi_r2
  alpha_pi <- conf$prior$alpha_pi
  a_phi <- conf$prior$a_phi
  alpha_delta <- conf$prior$alpha_delta
  alpha_zeta <- conf$prior$alpha_zeta

  p <- list()
  pi <- rdirichlet(J, alpha_pi)
  p$pi <- pi
  p$kappa <- matrix(
    data = c(qnorm(pi[, 1]), qnorm(pi[, 2] + pi[, 1])),
    nrow = J, ncol = 2
  )
  p$r2 <- rbeta(
    n = 1,
    shape1 = mu_r2 * phi_r2,
    shape2 = (1 - mu_r2) * phi_r2
  )
  p$t2 <- p$r2 / (1 - p$r2)
  p$phi <- rdirichlet(1, rep(a_phi, K * 4))
  p$beta <- rnorm(K, 0, p$phi[1:K] * p$t2)
  p$lambda <- rnorm(K, 0, p$phi[(K + 1):(2 * K)] * p$t2)
  p$gamma <- p$omega <- matrix(NA,
    nrow = J,
    ncol = K
  )
  p$delta <- p$zeta <- array(dim = c(J, K, 2))
  for (j in 1:J) {
    p$gamma[j, ] <- rnorm(K, 0, p$phi[((2 * K) + 1):(3 * K)] * p$t2)
    p$omega[j, ] <- rnorm(K, 0, p$phi[((3 * K) + 1):(4 * K)] * p$t2)
    p$delta[j, , ] <- rdirichlet(K, alpha_delta)
    p$zeta[j, , ] <- rdirichlet(K, alpha_zeta)
  }
  return(p)
}

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
