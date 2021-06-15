// Using both latent traits to predict a binary outcome (probit)
functions {
  int sum1d(int[] a, int x) {
    int s = 0;
    for (i in 1:num_elements(a)) {
      if (a[i] == x) s += 1;
    }
    return s;
  }
}
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
  int<lower = -1, upper = 1> outcome[Np];
}
transformed data {
  int n1 = sum1d(outcome, 1);
  int n0 = sum1d(outcome, 0);
  int nmiss = sum1d(outcome, -1);
  int idx1[n1];
  int idx0[n0];
  int idx_miss[nmiss];

  {
    int i1 = 0;
    int i0 = 0;
    int imiss = 0;
    for (i in 1:Np) {
      if (outcome[i] == 1) {
        i1 += 1;
        idx1[i1] = i;
      } else if (outcome[i] == 0) {
        i0 += 1;
        idx0[i0] = i;
      } else {
        imiss += 1;
        idx_miss[imiss] = i;
      }
    }
  }
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
  cholesky_factor_corr[3] R_chol;
  matrix[Np, 2] Eta_pr;
  real outcome_int;
  vector<lower = 0>[n1] z1;
  vector<upper = 0>[n0] z0;
  vector[nmiss] z_miss;
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
  outcome_int ~ normal(0, 2.5);

  lambda_p ~ lognormal(lambda_median, lambda_scale);
  sigma_eta ~ lognormal(lambda_median, lambda_scale);

  R_chol ~ lkj_corr_cholesky(shape_r);

  {
    vector[NpNi] p_s;
    vector[NpNi] disp_s;
    vector[NpNi] bb_a;
    vector[NpNi] bb_b;
    vector[NpNi - count_1 - count_max] lp_upper;
    vector[NpNi - count_1 - count_max] lp_lower;
    real lp_1;
    real lp_max;
    vector[3] Eta[Np];

    Eta[:, 1] = to_array_1d(Eta_pr[:, 1]);
    Eta[:, 2] = to_array_1d(Eta_pr[:, 2]);
    Eta[idx1, 3] = to_array_1d(z1 + outcome_int);
    Eta[idx0, 3] = to_array_1d(z0 + outcome_int);
    Eta[idx_miss, 3] = to_array_1d(z_miss + outcome_int);

    Eta ~ multi_normal_cholesky(rep_vector(0, 3), R_chol);

    p_s = inv_logit(p_lgt[i_ids] + lambda[i_ids] .* to_vector(Eta[p_ids, 1]));
    disp_s = exp(-r_lgt[i_ids] - sigma_eta * to_vector(Eta[p_ids, 2]));
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
  matrix[3, 3] R = multiply_lower_tri_self_transpose(R_chol);
  vector[Ni] p = inv_logit(p_lgt);
  vector[Ni] rho = inv_logit(r_lgt);
  real prop = Phi(-outcome_int);
  vector[Ni] scale = sqrt(N * p .* (1 - p) .* (1 + (N - 1) * rho));
  vector[2] coefs = mdivide_left_spd(R[1:2, 1:2], R[1:2, 3]);
}
