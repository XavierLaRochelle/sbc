library(extraDistr)
set.seed(1)
iter <- 200
N <- 1000
K <- 23
J <- 15
alpha_phi <- rep(1, K * 4)
precision <- 2
alpha_zeta <- alpha_delta <- precision / 2
p_x <- p_w <- c(1 / 3, 1 / 3, 1 / 3)

res <- c()
mu <- c()
X <- matrix(NA, N, K)
W <- matrix(NA, N, K)
beta <- matrix(NA, J, K)
lambda <- matrix(NA, J, K)
c <- sqrt(12 / 5)
for (i in 1:iter) {
  delta <- rbeta(K, alpha_delta, alpha_delta)
  zeta <- rbeta(K, alpha_zeta, alpha_zeta)
  jj <- sample(1:J, N, replace = T)
  for (k in 1:K) {
    X[, k] <- sample(c(0, delta[k], 1) * c, N, replace = T, prob = p_x)
    W[, k] <- sample(c(0, zeta[k], 1) * c, N, replace = T, prob = p_w)
  }

  phi <- rdirichlet(1, alpha_phi)
  Beta <- rnorm(K, 0, sqrt(phi[1:K]))
  Lambda <- rnorm(K, 0, sqrt(phi[(K + 1):(2 * K)]))
  for (k in 1:K) {
    beta[, k] <- rnorm(J, 0, sqrt(phi[(2 * K + k)]))
    lambda[, k] <- rnorm(J, 0, sqrt(phi[(3 * K + k)]))
  }
  for (n in 1:N) {
    mu[n] <- X[n, ] %*% (Beta + beta[jj[n], ]) +
      (X[n, ] * W[n, ]) %*% (Lambda + lambda[jj[n], ])
  }
  res <- c(res, mu)
}


var(res)
plot(density(res))
curve(dnorm(x), from = -4, to = 4, add = T, lty = 2)
