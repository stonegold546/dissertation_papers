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
}
parameters {
  real m_p_lgt;
  real<lower = 0> sd_p_lgt;
  vector[Ni] p_lgt;
  vector<lower = 0>[Ni] lambda_p;
  vector[Ni] r_lgt;
  vector[Np] eta;
}
transformed parameters {
  vector[Ni] lambda = lambda_p;

  lambda[neg_l] *= -1.0;
}
model {
  m_p_lgt ~ normal(0, 2);
  sd_p_lgt ~ normal(0, 2);
  p_lgt ~ normal(m_p_lgt, sd_p_lgt);
  r_lgt ~ normal(0, 2.5);

  lambda_p ~ lognormal(lambda_median, lambda_scale);
  eta ~ std_normal();

  {
    vector[NpNi] p_s = inv_logit(p_lgt[i_ids] + lambda[i_ids] .* eta[p_ids]);
    vector[NpNi] disp_s = exp(-r_lgt)[i_ids];
    n_1d ~ beta_binomial(N, p_s .* disp_s, (1 - p_s) .* disp_s);
  }
}
generated quantities {
  vector[Ni] p = inv_logit(p_lgt);
  vector[Ni] rho = inv_logit(r_lgt);
  vector[Ni] scale = sqrt(N * p .* (1 - p) .* (1 + (N - 1) * rho));
}
