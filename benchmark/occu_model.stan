// Standard single-species occupancy model (marginalized over latent z)
// Used for benchmarking against INLAocc and spOccupancy

data {
  int<lower=1> N;                          // number of sites
  int<lower=1> J;                          // visits per site
  array[N, J] int<lower=0, upper=1> y;     // detection history
  int<lower=1> K_occ;                      // occupancy covariates (incl. intercept)
  int<lower=1> K_det;                      // detection covariates (incl. intercept)
  matrix[N, K_occ] X_occ;                 // occupancy design matrix
  matrix[N * J, K_det] X_det;             // detection design matrix (stacked row-major)
}

parameters {
  vector[K_occ] beta_occ;
  vector[K_det] beta_det;
}

model {
  beta_occ ~ normal(0, 2.5);
  beta_det ~ normal(0, 2.5);

  for (i in 1:N) {
    real logit_psi = dot_product(X_occ[i], beta_occ);
    int base = (i - 1) * J;

    if (sum(y[i]) > 0) {
      // Detected at least once: z[i] = 1 with certainty
      target += log_inv_logit(logit_psi);
      for (j in 1:J) {
        target += bernoulli_logit_lpmf(y[i, j] | dot_product(X_det[base + j], beta_det));
      }
    } else {
      // Never detected: marginalize over z in {0, 1}
      real lp_present = log_inv_logit(logit_psi);
      for (j in 1:J) {
        lp_present += log1m_inv_logit(dot_product(X_det[base + j], beta_det));
      }
      target += log_sum_exp(lp_present, log1m_inv_logit(logit_psi));
    }
  }
}

generated quantities {
  vector[N] psi;
  for (i in 1:N) {
    psi[i] = inv_logit(dot_product(X_occ[i], beta_occ));
  }
}
