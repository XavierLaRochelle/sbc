functions {
  real induced_dirichlet_lpdf(vector kappa, vector alpha, real phi) {
    int K = num_elements(kappa) + 1;
    vector[K - 1] sigma = inv_logit(phi - kappa);
    vector[K] p;
    matrix[K, K] J = rep_matrix(0, K, K);

    // Induced ordinal probabilities
    p[1] = 1 - sigma[1];
    for (k in 2:(K - 1))
      p[k] = sigma[k - 1] - sigma[k];
    p[K] = sigma[K - 1];

    // Baseline column of Jacobian
    for (k in 1:K) J[k, 1] = 1;

    // Diagonal entries of Jacobian
    for (k in 2:K) {
      real rho = sigma[k - 1] * (1 - sigma[k - 1]);
      J[k, k] = - rho;
      J[k - 1, k] = rho;
    }

    return   dirichlet_lpdf(p | alpha)
      + log_determinant(J);
  }
  matrix mo(array[,] vector scale, array[,] int X, array[] int j) {

    int K = size(scale);
    int N = size(X[,1]);

    matrix[N, K] eta;
    for (n in 1:N) {
      for (k in 1:K) {
        if (X[n, k] == 1) {
          eta[n, k] = 0;
        } else {
          eta[n, k] = sum(scale[j[n], k, 1:(X[n, k]) - 1]);
        }
      }
    }
    return(eta);
  }
}

data {
  int N;               // Number of observations (patients)
  int J;               // Number of levels (evaluators)
  int K;               // Number of predictors
  int D_y;                   // Number of categories of the outcome
  int D_X;                   // Number of categories of the predictors
  array[N] int y;     // Observed ordinal outcome (risk estimates)
  array[N] int j;     // Level (evaluator) index
  array[N, K] int X;  // predictor matrix
}
transformed data {
  real<lower = 0> half_K = 0.5 * K;
  real<lower = 0> sqrt_Nm1 = sqrt(N - 1.0);
}
parameters {

  array[J] simplex[D_y] pi;             // Probabilities of the data
  array[J] unit_vector[K] u;            // Latent effects of the predictors per level
  array[J,K] simplex[D_X - 1] delta;    // Normalized distances across categories of the predictors per level
}
transformed parameters {
  array[J] vector[K] beta;
  array[J] vector[D_y - 1] cutpoints;
  {
    real delta_y = inv_sqrt(1 - R2);
    beta = u * sqrt(R2) * delta_y * sqrt_Nm1;
    cutpoints = logit(cumulative_sum(pi[1:(D_y - 1)]) * delta_y);
  }
}
model {

  // Prior
  for(n in 1:N) {
    phi[n] = mo(X) * beta[j[n]] ;
  }

  for(n in 1:N) {
    y[n] ~ ordered_logistic(phi[n], kappa[j[n]]);
  }
}

