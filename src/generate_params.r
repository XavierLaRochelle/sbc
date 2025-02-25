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
