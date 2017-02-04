functions {
  #include 'common_functions.stan'
}

data {
  int<lower=1> N; // number of bins without padding
  int<lower=4> Ntot; // number of bins with padding
  int<lower=1> N_vis; // number of bins we observe for each process
  int<lower=1, upper=Ntot/2 - 1> n_nonzero_freq; // number of non-zero frequencies to use
  int<lower=0> counts[N];
  vector[N] offset;
  real<lower=0> delta;
  int<lower=0> npars;
  int<lower=0, upper=1> GMO_FLAG;
  vector[npars * GMO_FLAG] fixed_phi;
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
  lengthscale ~ normal(0, 10);
  target += log(lengthscale);
  sigma ~ lognormal(0, 1);
  target += log(sigma);
  mu ~ normal(0, 5);
  nu ~ normal(0, 1);
  target += log(nu);

  counts[1:N_vis] ~ poisson_log(x[1:N_vis] * sigma + mu + log(offset[1:N_vis]));
  
}

generated quantities {
  vector[N] rates;
  rates = exp(x[1:N] * sigma + mu);
}
