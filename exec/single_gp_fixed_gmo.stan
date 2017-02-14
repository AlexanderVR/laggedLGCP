functions {
  #include 'common_functions.stan'
}

data {
  int<lower=1> N; // number of bins without padding
  int<lower=4> Ntot; // number of bins with padding
  int<lower=1> N_vis; // number of bins we have observed data for
  int<lower=1, upper=Ntot/2 - 1> n_nonzero_freq; // number of non-zero frequencies to use
  int<lower=-1> counts[N];
  vector[N] offset;
  int<lower=0, upper=1> coalescent;
  real<lower=0> prior_lengthscale_mean;
  real<lower=0> prior_smoothness_mean;
  real<lower=0> prior_smoothness_std;
  real<lower=0> delta;
  int<lower=0> npars;
  int<lower=0, upper=1> GMO_FLAG;
  vector[npars * GMO_FLAG] fixed_phi;
}

transformed data {
  int idx[N_vis];  // indicies of observed data
  int loc;
  real signs;
  // coalescence rate is proportional to 1/N_e(t)
  signs = (coalescent==1) ? -1 : 1;
  loc = 1;
  for (k in 1:N) {
    if (counts[k] >= 0) {
      idx[loc] = k;
      loc = loc + 1;
    }
  }
}

parameters {
  vector[npars * (1 - GMO_FLAG)] phi;
  vector[n_nonzero_freq] im;
  vector[n_nonzero_freq] re;
}

transformed parameters {
  vector[Ntot] x; // latent correlated GP
  real sigma;
  real mu;
  real lengthscale;
  real nu;
  
  if (GMO_FLAG) {
    mu = fixed_phi[1];
    sigma = exp(fixed_phi[2]);
    lengthscale = exp(fixed_phi[3]);
    nu = exp(fixed_phi[4]);
  }
  else {
    mu = phi[1];
    sigma = exp(phi[2]);
    lengthscale = exp(phi[3]);
    nu = exp(phi[4]);
  }
  x = transform_to_matern(re, im, Ntot, nu, lengthscale, delta, 0.);
}

model {
  to_vector(re) ~ normal(0, .5);
  to_vector(im) ~ normal(0, .5);
  lengthscale ~ lognormal(log(prior_lengthscale_mean), 1);
  target += log(lengthscale);
  sigma ~ lognormal(0, 1);
  target += log(sigma);
  mu ~ normal(0, 5);
  nu ~ normal(prior_smoothness_mean, prior_smoothness_std);
  target += log(nu);
  
  counts[idx] ~ poisson_log(signs * (x[idx] * sigma + mu) + log(offset[idx]));
  
}

generated quantities {
  vector[N] rates;
  rates = exp(x[1:N] * sigma + mu);
}
