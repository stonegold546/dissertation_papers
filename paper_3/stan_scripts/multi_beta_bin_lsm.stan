data {
  int Np;
  int Ni;
  int NpNi;
  int n_1d[NpNi];
  int p_ids[NpNi];
  int i_ids[NpNi];
  int N;
  real lambda_median;
  real<lower = 0> lambda_scale;
  int n_neg;
  int<lower = 1, upper = Ni> neg_l[n_neg];
  real<lower = 0> shape_r;
}
parameters {
  real m_p_lgt;
  real<lower = 0> sd_p_lgt;
  vector[Ni] p_lgt;
  vector<lower = 0>[Ni] lambda_p;
  real m_r_lgt;
  real<lower = 0> sd_r_lgt;
  vector[Ni] r_lgt;
  real<lower = 0> sigma_eta;
  cholesky_factor_corr[2] R_chol;
  vector[2] Eta[Np];
}
transformed parameters {
  vector[Ni] lambda = lambda_p;

  lambda[neg_l] *= -1.0;
}
model {
  m_p_lgt ~ normal(0, 2);
  sd_p_lgt ~ normal(0, 2);
  p_lgt ~ normal(m_p_lgt, sd_p_lgt);
  m_r_lgt ~ normal(0, 2);
  sd_r_lgt ~ normal(0, 2);
  r_lgt ~ normal(m_r_lgt, sd_r_lgt);

  lambda_p ~ lognormal(lambda_median, lambda_scale);
  sigma_eta ~ student_t(3, 0, 1);

  R_chol ~ lkj_corr_cholesky(shape_r);
  Eta ~ multi_normal_cholesky(rep_vector(0, 2), R_chol);

  {
    vector[NpNi] p_s = inv_logit(p_lgt[i_ids] + lambda[i_ids] .* to_vector(Eta[p_ids, 1]));
    vector[NpNi] disp_s = exp(-r_lgt[i_ids] - sigma_eta * to_vector(Eta[p_ids, 2]));
    n_1d ~ beta_binomial(N, p_s .* disp_s, (1 - p_s) .* disp_s);
  }
}
generated quantities {
  real r = R_chol[2, 1];
  vector[Ni] p = inv_logit(p_lgt);
  vector[Ni] rho = inv_logit(r_lgt);
  vector[Ni] scale = sqrt(N * p .* (1 - p) .* (1 + (N - 1) * rho));
}
