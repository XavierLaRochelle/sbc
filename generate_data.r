library(extraDistr)

# Step 1: choose configuration
conf <- list(
  seed = 10,
  N = 100,
  J = 5
)

# Set 2: Choose the prior
conf$prior <- list(
  mu_r2 = 0.1,
  phi_r2 = 0.5,
  a_phi = 1,
  alpha_pi = c(1, 1, 1),
  alpha_delta = c(1, 1),
  alpha_zeta = c(1, 1)
)

# Step 3: Generate the parameters from the prior
generate_params <- function(conf) {
  set.seed(conf$seed)
  J <- conf$J
  mu_r2 <- conf$prior$mu_r2
  phi_r2 <- conf$prior$phi_r2
  alpha_pi <- conf$prior$alpha_pi
  a_phi <- conf$prior$a_phi
  alpha_delta <- conf$prior$alpha_delta
  alpha_zeta <- conf$prior$alpha_zeta

  p <- list()
  p$pi <- rdirichlet(J, alpha_pi)
  p$kappa <- matrix(
    data = c(qnorm(p$pi[, 1]), qnorm(p$pi[, 2] + p$pi[, 1])),
    nrow = J, ncol = 2
  )
  p$r2 <- rbeta(
    n = 1,
    shape1 = mu_r2 * phi_r2,
    shape2 = (1 - mu_r2) * phi_r2
  )
  p$t2 <- p$r2 / (1 - p$r2)
  p$phi <- rdirichlet(1, rep(a_phi, 23 * 4))
  p$beta <- rnorm(23, 0, p$phi[1:23] * p$t2)
  p$lambda <- rnorm(23, 0, p$phi[24:46] * p$t2)
  p$gamma <- p$omega <- matrix(NA, nrow = J, ncol = 23)
  p$delta <- p$zeta <- array(dim = c(J, 23, 2))
  for (j in 1:J) {
    p$gamma[j, ] <- rnorm(23, 0, p$phi[47:69] * p$t2)
    p$omega[j, ] <- rnorm(23, 0, p$phi[70:92] * p$t2)
    p$delta[j, , ] <- rdirichlet(23, alpha_delta)
    p$zeta[j, , ] <- rdirichlet(23, alpha_zeta)
  }
  return(p)
}

params <- generate_params(conf)
# Step 4: Generate the data
generate_data <- function(conf, params) {
  set.seed(conf$seed)
  N <- conf$N
  J <- conf$J


  d <- list()
  return(d)
}

generate_data(conf, params)
