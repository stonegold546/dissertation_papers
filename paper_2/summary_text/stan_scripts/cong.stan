// congeneric model (heteroskedastic) with log-normal loadings
data {
  int<lower = 2> Np;
  int<lower = 1> Ni;
  int<lower = 2> N;
  real<lower = 0> alpha_scale;
  real<lower = 0> beta_scale;
  real lambda_median;
  real<lower = 0> lambda_scale;
  real<lower = 0> sigma_scale;
  int<lower = 1, upper = Np> resp_id[N];
  int<lower = 1, upper = Ni> item_id[N];
  vector[N] y;
  int ret_yhat;
  int ret_ll;
}
transformed data {
  int Ny = 0;
  int Nll = 0;

  if (ret_yhat == 1) Ny = N;
  if (ret_ll == 1) Nll = N;
}
parameters {
  real alpha;
  vector[Ni] beta_dev;
  vector<lower = 0>[Ni] lambda;
  vector[Np] theta_p;
  vector<lower = 0>[Ni] sigma;
}
model {
  alpha ~ normal(0, alpha_scale);
  beta_dev ~ normal(0, beta_scale);

  lambda ~ lognormal(lambda_median, lambda_scale);
  theta_p ~ std_normal();

  sigma ~ normal(0, sigma_scale);

  {
    vector[N] sigma_y;
    vector[N] mu;
    for (i in 1:N) {
      sigma_y[i] = sigma[item_id[i]];
      mu[i] = alpha + beta_dev[item_id[i]] + lambda[item_id[i]] * theta_p[resp_id[i]];
    }
    y ~ normal(mu, sigma_y);
  }
}
generated quantities {
  vector[Ni] beta = alpha + beta_dev;
  vector[Nll] log_lik;
  vector[Ny] yhat;

  {
    real mu;
    for (i in 1:max(Nll, Ny)) {
      mu = alpha + beta_dev[item_id[i]] + lambda[item_id[i]] * theta_p[resp_id[i]];
      if (Nll > 0) log_lik[i] = normal_lpdf(y[i] | mu, sigma[item_id[i]]);
      if (Ny > 0) yhat[i] = normal_rng(mu, sigma[item_id[i]]);
    }
  }
}
