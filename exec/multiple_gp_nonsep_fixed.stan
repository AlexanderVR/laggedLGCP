functions {
	#include 'common_functions.stan'
}

data {
  int<lower=1> N; // number of bins without padding
  int<lower=4> Ntot; // number of bins to use
  int<lower=1> n_series; // number of processes
  int<lower=1> N_vis; // total number of observed bins
  int<lower=1, upper=Ntot/2 - 1> n_nonzero_freq; // number of non-zero frequencies to use
  int<lower=-1> counts[N, n_series];
  matrix[N, n_series] offset;
  real<lower=0> max_lag;
  real<lower=0> delta;
  int<lower=0, upper=1> coalescent;
  vector<lower=0>[n_series] lengthscale;
  vector<lower=0>[n_series] nu;
  vector[n_series] mu;
  vector<lower=0>[n_series] sigma;
}

transformed data {
  int loc;
  int pos[n_series + 1]; // position of each group
  int idx[N_vis];  // indicies of observed data under ragged array structure

  vector[n_series] signs;
  for (j in 1:n_series)
    signs[j] = (j==1 && coalescent==1) ? -1 : 1;
  
  // find indicies for non-negative elements of "counts" input 
  // active indicies for the j'th series are:
  // idx[pos[j] : pos[j + 1] - 1].
  loc = 1;
  pos[1] = 1;
  for (j in 1:n_series) {
    for (k in 1:N) {
      if (counts[k, j] >= 0) {
        idx[loc] = k;
        loc = loc + 1;  
      }
    }
    pos[j + 1] = loc;
  }
}

parameters {
  cholesky_factor_corr[n_series] cor;
  vector<lower=-max_lag, upper=max_lag>[n_series - 1] shifts;
  matrix[n_nonzero_freq, n_series] re;
  matrix[n_nonzero_freq, n_series] im;
}

transformed parameters {
  matrix[Ntot, n_series] x; // latent correlated GP

  for (j in 1:n_series) {
    real shift;

    shift = (j == n_series) ? 0. : shifts[j]; // last process has shift = 0.
    x[, j] = transform_to_matern(re[, j], im[, j], Ntot, nu[j], lengthscale[j], delta, shift);
  }

}
model {
  for (n in 1:n_nonzero_freq) {
    re[n, ] ~ multi_normal_cholesky(rep_vector(0, n_series), cor / sqrt(2));
    im[n, ] ~ multi_normal_cholesky(rep_vector(0, n_series), cor / sqrt(2));
  }
  
  cor ~ lkj_corr_cholesky(1.);
  
  for (j in 1:n_series) {
    // coalescence rate is proportional to 1/N_e(t)
    // This would be much cleaner if Stan allowed ragged array structures
    counts[idx[pos[j] : (pos[j + 1] - 1)], j] ~ poisson_log(signs[j] * (x[idx[pos[j] : (pos[j + 1] - 1)], j] * sigma[j] + mu[j]) + log(offset[idx[pos[j] : (pos[j + 1] - 1)], j]));
  }
}

generated quantities {
  matrix[N, n_series] rates;
  matrix[n_series, n_series] correlations;

  rates = exp(diag_post_multiply(x[1:N, ], sigma) + rep_matrix(mu', N));
  correlations = cor * cor';
}
