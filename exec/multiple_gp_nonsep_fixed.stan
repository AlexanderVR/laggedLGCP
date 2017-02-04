functions {
	#include 'common_functions.stan'
}

data {
  int<lower=1> N; // number of bins without padding
  int<lower=4> Ntot; // number of bins to use
  int<lower=1> n_series; // number of processes
  int<lower=1> N_vis[n_series]; // number of bins we observe for each process
  int<lower=1, upper=Ntot/2 - 1> n_nonzero_freq; // number of non-zero frequencies to use
  int<lower=0> counts[N, n_series];
  matrix[N, n_series] offset;
  real<lower=0> max_lag;
  real<lower=0> delta;
  int<lower=0, upper=1> coalescent;
  vector<lower=0>[n_series] lengthscale;
  vector<lower=0>[n_series] nu;
  vector[n_series] mu;
  vector<lower=0>[n_series] sigma;
}

parameters {
  cholesky_factor_corr[n_series] cor;
  vector[n_series - 1] shifts;
  matrix[n_nonzero_freq, n_series] re;
  matrix[n_nonzero_freq, n_series] im;
}

transformed parameters {
  matrix[Ntot, n_series] x; // latent correlated GP

  for (j in 1:n_series) {
    real shift;
    real signs;
    if (coalescent && (j == 1))
      signs = -1; // coalescence rate is proportional to 1/N_e(t)
    else
      signs = 1;  

    shift = (j == n_series) ? 0. : shifts[j]; // last process has shift = 0.
    x[, j] = signs * transform_to_matern(re[, j], im[, j], Ntot, nu[j], lengthscale[j], delta, shift);
  }

}
model {
  for (n in 1:n_nonzero_freq) {
    re[n, ] ~ multi_normal_cholesky(rep_vector(0, n_series), cor / sqrt(2));
    im[n, ] ~ multi_normal_cholesky(rep_vector(0, n_series), cor / sqrt(2));
  }
  
  shifts ~ normal(0, max_lag / 2.);
  cor ~ lkj_corr_cholesky(1.);
  
  for (j in 1:n_series) {
    counts[1:N_vis[j], j] ~ poisson_log(x[1:N_vis[j], j] * sigma[j] + mu[j] + log(offset[1:N_vis[j], j]));
  }
}

generated quantities {
  matrix[N, n_series] rates;
  matrix[n_series, n_series] correlations;

  rates = exp(diag_post_multiply(x[1:N, ], sigma) + rep_matrix(mu', N));
  correlations = cor * cor';
}
