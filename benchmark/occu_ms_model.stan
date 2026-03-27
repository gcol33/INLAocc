// Multi-species occupancy model with community hyperpriors
// Marginalized over latent z (HMC can't sample discrete states)

data {
  int<lower=1> S;                              // species
  int<lower=1> N;                              // sites
  int<lower=1> J;                              // visits per site
  array[S, N, J] int<lower=0, upper=1> y;      // detection data
  int<lower=1> K_occ;                          // occupancy covariates (incl. intercept)
  int<lower=1> K_det;                          // detection covariates (incl. intercept)
  matrix[N, K_occ] X_occ;                     // occupancy design matrix
  matrix[N * J, K_det] X_det;                 // detection design matrix (stacked)
}

parameters {
  vector[K_occ] mu_beta;                       // community occupancy means
  vector[K_det] mu_alpha;                      // community detection means
  vector<lower=0>[K_occ] sigma_beta;           // community occupancy SDs
  vector<lower=0>[K_det] sigma_alpha;          // community detection SDs
  array[S] vector[K_occ] beta;                // species occupancy coefficients
  array[S] vector[K_det] alpha;               // species detection coefficients
}

model {
  mu_beta ~ normal(0, 2.5);
  mu_alpha ~ normal(0, 2.5);
  sigma_beta ~ exponential(1);
  sigma_alpha ~ exponential(1);

  for (s in 1:S) {
    beta[s] ~ normal(mu_beta, sigma_beta);
    alpha[s] ~ normal(mu_alpha, sigma_alpha);

    for (i in 1:N) {
      real logit_psi = dot_product(X_occ[i], beta[s]);
      int base = (i - 1) * J;

      if (sum(y[s, i]) > 0) {
        target += log_inv_logit(logit_psi);
        for (j in 1:J) {
          target += bernoulli_logit_lpmf(y[s, i, j] |
                      dot_product(X_det[base + j], alpha[s]));
        }
      } else {
        real lp = log_inv_logit(logit_psi);
        for (j in 1:J) {
          lp += log1m_inv_logit(dot_product(X_det[base + j], alpha[s]));
        }
        target += log_sum_exp(lp, log1m_inv_logit(logit_psi));
      }
    }
  }
}

generated quantities {
  array[S] vector[N] psi;
  for (s in 1:S) {
    for (i in 1:N) {
      psi[s][i] = inv_logit(dot_product(X_occ[i], beta[s]));
    }
  }
}
