source("./src/utils.r")

# Step 1: choose a configuration
conf <- list(
  seed = 10,
  N = 1000,
  J = 5,
  K = 23
)

# Set 2: Choose the prior
prior <- list(
  mu_r2 = 0.1,
  phi_r2 = 0.5,
  a_phi = 1,
  alpha_pi = c(1, 1, 1),
  alpha_delta = c(1, 1),
  alpha_zeta = c(1, 1)
)

# Step 3: Generate parameters from the prior
params <- generate_params(conf, prior)

# Step 4: Generate data from the parameters
data <- generate_data(conf, params)

# Step 5: Generate posterior draws from the data
draws <- generate_draws(conf, data)

obj <- list(conf = conf, params = params, data = data)
saveRDS(object = obj, file = "./data/temp/test.rds")
