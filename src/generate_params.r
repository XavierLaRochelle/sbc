generate_params <- function(seed, J, K, mu_r2, phi_r2, alpha_pi, a_phi, alpha_delta, alpha_zeta) {
  require(extraDistr)
  set.seed(seed)

  p <- list()
  pi <- rdirichlet(J, alpha_pi)
  p$kappa <- matrix(
    data = c(qnorm(pi[, 1]), qnorm(pi[, 2] + pi[, 1])),
    nrow = J, ncol = 2
  )
  r2 <- rbeta(
    n = 1,
    shape1 = mu_r2 * phi_r2,
    shape2 = (1 - mu_r2) * phi_r2
  )
  t2 <- r2 / (1 - r2)
  phi <- rdirichlet(1, rep(a_phi, K * 4))
  p$beta <- rnorm(K, 0, phi[1:K] * t2)
  p$lambda <- rnorm(K, 0, phi[(K + 1):(2 * K)] * t2)
  p$gamma <- p$omega <- matrix(NA,
    nrow = J,
    ncol = K
  )
  p$delta <- p$zeta <- array(dim = c(J, K, 2))
  for (j in 1:J) {
    p$gamma[j, ] <- rnorm(K, 0, phi[((2 * K) + 1):(3 * K)] * t2)
    p$omega[j, ] <- rnorm(K, 0, phi[((3 * K) + 1):(4 * K)] * t2)
    p$delta[j, , ] <- rdirichlet(K, alpha_delta)
    p$zeta[j, , ] <- rdirichlet(K, alpha_zeta)
  }
  return(p)
}
