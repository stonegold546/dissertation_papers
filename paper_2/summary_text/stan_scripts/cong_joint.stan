// congeneric model (heteroskedastic) with log-normal loadings
data {
  int<lower = 2> Np;
  int<lower = 1> Ni;
  int<lower = 2> N;
  real<lower = 0> alpha_scale;
  real<lower = 0> beta_scale;
  real lambda_median;
  real<lower = 0> lambda_scale;
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
  real ln_alpha;
  vector[Ni] ln_sigma2_dev;
  real<lower = 0> ls_scale;
  real<lower = 0> ll_scale;
  vector<lower = 0>[2] ln_lambda_pos;
  vector[Ni - 2] ln_lambda_oth;
  vector[2] Theta[Np];
  real<lower = 0, upper = 1> r_tr;
}
transformed parameters {
  // first and third elements assumed positive based on figure 5 in paper
  vector[Ni] ln_lambda;
  ln_lambda[1] = ln_lambda_pos[1];
  ln_lambda[3] = ln_lambda_pos[2];
  ln_lambda[2] = ln_lambda_oth[1];
  ln_lambda[4] = ln_lambda_oth[2];
  ln_lambda[5] = ln_lambda_oth[3];
  ln_lambda[6] = ln_lambda_oth[4];
}
model {
  matrix[2, 2] R;

  R[1, 1] = 1;
  R[2, 2] = 1;
  R[1, 2] = r_tr * 2 - 1;
  R[2, 1] = r_tr * 2 - 1;

  r_tr ~ beta(2, 2);

  alpha ~ normal(0, alpha_scale);
  beta_dev ~ normal(0, beta_scale);

  lambda ~ lognormal(lambda_median, lambda_scale);

  // prior predictive checking would yield informative priors here
  ln_alpha ~ normal(0, 5);  // vague prior, do not limit scale of overall variance mean
  ln_sigma2_dev ~ normal(0, ls_scale);  // using hyperprior
  ls_scale ~ std_normal();

  ln_lambda ~ normal(0, ll_scale);  // using hyperprior
  ll_scale ~ std_normal();

  Theta ~ multi_normal(rep_vector(0, 2), R);

  {
    vector[N] sigma_y;
    vector[N] mu;
    for (i in 1:N) {
      sigma_y[i] = sqrt(exp(ln_alpha + ln_sigma2_dev[item_id[i]] + ln_lambda[item_id[i]] * Theta[resp_id[i], 2]));
      mu[i] = alpha + beta_dev[item_id[i]] + lambda[item_id[i]] * Theta[resp_id[i], 1];
    }
    y ~ normal(mu, sigma_y);
  }
}
generated quantities {
  vector[Ni] beta = alpha + beta_dev;
  vector[Nll] log_lik;
  vector[Ny] yhat;
  vector[Ni] sigma = sqrt(exp(ln_alpha + ln_sigma2_dev));
  real r = r_tr * 2 - 1;

  {
    real sigma_y;
    real mu;
    for (i in 1:max(Nll, Ny)) {
      sigma_y = sqrt(exp(ln_alpha + ln_sigma2_dev[item_id[i]] + ln_lambda[item_id[i]] * Theta[resp_id[i], 2]));
      mu = alpha + beta_dev[item_id[i]] + lambda[item_id[i]] * Theta[resp_id[i], 1];
      if (Nll > 0) log_lik[i] = normal_lpdf(y[i] | mu, sigma_y);
      if (Ny > 0) yhat[i] = normal_rng(mu, sigma_y);
    }
  }
}
