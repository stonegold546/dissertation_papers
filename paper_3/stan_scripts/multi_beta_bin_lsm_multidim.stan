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
  real<lower = 0> shape_r;
  int Nl;
  int Nf;
  int Load_Pattern[Ni, Nf];
}
transformed data {
  int Nfp1 = Nf + 1;
}
parameters {
  real m_p_lgt;
  real<lower = 0> sd_p_lgt;
  vector[Ni] p_lgt;
  vector<lower = 0>[Nl] lambda_p;
  real m_r_lgt;
  real<lower = 0> sd_r_lgt;
  vector[Ni] r_lgt;
  real<lower = 0> sigma_eta;
  cholesky_factor_corr[Nfp1] R_chol;
  vector[Nfp1] Eta[Np];
}
transformed parameters {
  matrix[Ni, Nf] Lambda_mat = rep_matrix(0, Ni, Nf);

  {
    int pos = 0;
    for (i in 1:Ni) {
      for (j in 1:Nf) {
        if (Load_Pattern[i, j] == 1) {
          pos += 1;
          Lambda_mat[i, j] = lambda_p[pos];
        } else if (Load_Pattern[i, j] == -1) {
          pos += 1;
          Lambda_mat[i, j] = lambda_p[pos] * -1;
        }
      }
    }
  }
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
  Eta ~ multi_normal_cholesky(rep_vector(0, Nfp1), R_chol);

  {
    vector[NpNi] p_s;
    vector[NpNi] disp_s;
    vector[NpNi] bb_a;
    vector[NpNi] bb_b;

    for (i in 1:NpNi) {
      p_s[i] = inv_logit(p_lgt[i_ids[i]] + Lambda_mat[i_ids[i], ] * Eta[p_ids[i], 1:Nf]);
    }
    disp_s = exp(-r_lgt[i_ids] - sigma_eta * to_vector(Eta[p_ids, Nfp1]));
    bb_a = p_s .* disp_s;
    bb_b = disp_s - bb_a;

    n_1d ~ beta_binomial(N, bb_a, bb_b);
  }
}
generated quantities {
  vector[Ni] p = inv_logit(p_lgt);
  vector[Ni] rho = inv_logit(r_lgt);
  vector[Ni] scale = sqrt(N * p .* (1 - p) .* (1 + (N - 1) * rho));
  matrix[Nfp1, Nfp1] R = multiply_lower_tri_self_transpose(R_chol);
}
