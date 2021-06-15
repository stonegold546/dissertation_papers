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
  int N_int;
  int cuts[N_int, 2];
  int count_1;
  int count_max;
  int pos_1[count_1];
  int pos_max[count_max];
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
  matrix[2, Np] Eta_pr;
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
  sigma_eta ~ lognormal(lambda_median, lambda_scale);

  R_chol ~ lkj_corr_cholesky(shape_r);
  to_vector(Eta_pr) ~ std_normal();

  {
    vector[NpNi] p_s;
    vector[NpNi] disp_s;
    vector[NpNi] bb_a;
    vector[NpNi] bb_b;
    vector[NpNi - count_1 - count_max] lp_upper;
    vector[NpNi - count_1 - count_max] lp_lower;
    real lp_1;
    real lp_max;
    matrix[Np, 2] Eta = (R_chol * Eta_pr)';

    p_s = inv_logit(p_lgt[i_ids] + lambda[i_ids] .* Eta[p_ids, 1]);
    disp_s = exp(-r_lgt[i_ids] - sigma_eta * Eta[p_ids, 2]);
    bb_a = p_s .* disp_s;
    bb_b = disp_s - bb_a;

    if (cuts[1, 2] == 0) lp_1 = beta_binomial_lpmf(0 | N, bb_a[pos_1], bb_b[pos_1]);
    else lp_1 = beta_binomial_lcdf(cuts[1, 2] | N, bb_a[pos_1], bb_b[pos_1]);

    lp_max = beta_binomial_lccdf(cuts[N_int, 1] - 1 | N, bb_a[pos_max], bb_b[pos_max]);

    {
      int pos = 0;
      for (i in 1:NpNi) {
        if (n_1d[i] > 1 && n_1d[i] < N_int) {
          pos += 1;
          lp_upper[pos] = beta_binomial_lcdf(cuts[n_1d[i], 2] | N, bb_a[i], bb_b[i]);
          lp_lower[pos] = beta_binomial_lcdf(cuts[n_1d[i], 1] - 1 | N, bb_a[i], bb_b[i]);
        }
      }
    }

    target += lp_1;
    target += lp_max;
    target += log_diff_exp(lp_upper, lp_lower);
  }
}
generated quantities {
  real r = R_chol[2, 1];
  vector[Ni] p = inv_logit(p_lgt);
  vector[Ni] rho = inv_logit(r_lgt);
  vector[Ni] scale = sqrt(N * p .* (1 - p) .* (1 + (N - 1) * rho));
}
