// R2D2M2 Version for simulation example


/* Version of R2D2M2 prior K grouping factors with multiple levels L 
   and D fixed effects and intercept.
* All groups have the same number of levels Lg as well as the same number of coefficients
* Dg. Dg is the number of group level effects, this is counting random intercepts as well.
* i.e Group k level l has Dg level effects.
* D is the number of fixed covariates and the intercept.
* This version is able to work with data where D, K, Lg can vary with no need to manually
* modify the Stan code. 
* We still need to handle the case when groups have different amount of levels and covariates 
* per group. 

* It is possible to briefly modify this code in order to work with
* q <= p random effects. The order of the random effects should consider the 
* relationship with the ith covariate, since we scale by sdx_i.

*/

functions {
  
  vector R2D2(vector z, vector phi, real tau2) {
    /* Efficient computation of the R2D2 prior
    * Args:
    *   z: standardized population-level coefficients
    *   phi: local weight parameters
    *   tau2: global scale parameter (sigma is inside tau2)
    * Returns:
      *   population-level coefficients following the R2D2 prior
    */
      return  z .* sqrt(phi * tau2);
  }

}
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  int<lower=1> D;  // number of population-level effects including intercept
  matrix[N, D] X;  // population-level design matrix including column of 1s
  int<lower=0> K; // number of groups
  
  //---- data for group-level effects 
  
  int<lower=1> Lg;  // number of  levels per group (constant)
  int<lower=1> Dg; // number of coefficients per level per group (D_g constant per group)
  int<lower=1> J[N,K]; // grouping indicator matrix per observation per group K
  
  //---- group-level predictor values 
  matrix[Dg,N] Z[K]; 
  
  
  //---- data for the R2D2 prior
  vector<lower=0>[ (D-1)+K+(Dg-1)*K] R2D2_alpha; 
  real<lower=0> R2D2_mean_R2;  // mean of the R2 prior
  real<lower=0> R2D2_prec_R2;  // precision of the R2 prior

}

transformed data {
  int Dc = D - 1; 
  matrix[N, Dc] Xc;  // centered version of X without an intercept
  vector[Dc] means_X;  // column means of X before centering
  vector[N] Yc;
  real Ymean;
  
  for (i in 2:D) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = (X[, i] - means_X[i - 1]) ; 
  }
  
  Ymean= mean(Y);
  for (i in 1:N) {
    Yc[i]= Y[i]-mean(Y);
  }
  
}

parameters {
  real Intercept;  // temporary intercept for centered predictors
  
  vector[Dc] zb; // standardized population-level effects
  matrix[Dg,Lg] z[K]; // standardized group-level effects 
  
  // local parameters for the R2D2M2 prior
  simplex[Dc+K+(Dg-1)*K] R2D2_phi; 
  // R2D2 shrinkage parameters 
  // Convention of indexing: First Dc for fixed effects, group of K for random intercepts,
  // Batches of Dc for each group.  
  real<lower=0,upper=1> R2D2_R2;  // R2 parameter
  real<lower=0> sigma;  // residual error
}

transformed parameters {
  
  vector[Dc] b;  // population-level effects
  matrix[Dg,Lg] r[K]; // actual group-level effects (includes random intercept)
  real R2D2_tau2;  // global R2D2 scale parameter
  
  R2D2_tau2 =  R2D2_R2 / (1 - R2D2_R2);
  
  // compute actual regression coefficients
  b = R2D2(zb, R2D2_phi[1:Dc], (sigma^2) * R2D2_tau2);  
  
  for(k in 1:K){
    // random intercepts
    // Dc+k is the kth random intercept
    r[k,1,] = (sigma * sqrt(R2D2_tau2 * R2D2_phi[Dc+k ]) * (z[k,1,]));  
    for(d in 2: Dg){
      // group level effects
      // (k-1)Dc indexes the beginning of the kth batch of scales
      //careful with sds_X
      r[k,d,]= sigma * sqrt(R2D2_tau2 * R2D2_phi[Dc+K+ (k-1)*(Dg-1) +(d-1) ]) * (z[k,d,]);  
    }
  }
}

model {
  // likelihood including constants
  
  
  // initialize linear predictor term
  vector[N] mu = Intercept + rep_vector(0.0, N);
  for (n in 1:N) {
    // add more terms to the linear predictor 
    for(k in 1:K){
      mu[n]+=dot_product(r[k,,J[n,k]], Z[k,,n]) ;
      }
  }
  target += normal_id_glm_lpdf(Yc | Xc, mu, b, sigma); // mu+ Xc*b

  // priors including constants
  target += beta_lpdf(R2D2_R2 | R2D2_mean_R2 * R2D2_prec_R2, (1 - R2D2_mean_R2) * R2D2_prec_R2); // R^2 
  target += dirichlet_lpdf(R2D2_phi | R2D2_alpha); //phi
  target += std_normal_lpdf(zb); //zb: overall effects
  
  for(k in 1:K){
    for(d in 1: Dg){
      target += std_normal_lpdf(z[k,d,]); // z: varying effects
    }
  }
  
  target += normal_lpdf(Intercept | 0, 5);  // Intercept
  target += student_t_lpdf(sigma | 3, 0, sd(Yc));  
  
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Ymean+Intercept - dot_product(means_X, b); 
  
}
