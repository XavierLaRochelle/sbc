source("./utils.r")
library(cmdstanr)

# Step 1: choose configuration
conf <- list(
  seed = 10,
  N = 1000,
  J = 5,
  K = 23
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
params <- generate_params(conf)

# Step 4: Generate the data
data <- generate_data(conf, params)

# Step 5: Generate posterior draws
draws <- generate_draws(conf, data)

obj <- list(conf = conf, params = params, data = data)
saveRDS(object = obj, file = "./data/test.rds")
