data {
  int Np;
  int Ni;
  int NpNi;
  int n[NpNi];
  int p_ids[NpNi];
  int i_ids[NpNi];
  real lambda_median;
  real<lower = 0> lambda_scale;
  int N_int;
  int cuts[N_int + 1];
  int count_1;
  int count_max;
  int pos_1[count_1];
  int pos_max[count_max];
  int Nl;
  int Nf;
  int Load_Pattern[Ni, Nf];
}
transformed data {
  int N = cuts[N_int + 1];
}
parameters {
  real mu_m;
  real<lower = 0> sigma_m;
  vector[Ni] means;
  vector<lower = 0>[Nl] lambda;
  cholesky_factor_corr[Nf] R_chol;
  vector[Nf] eta[Np];
  vector<lower = 0>[Ni] precisions;
  real<lower = 0> gamma_a;
  real<lower = 0> gamma_b;
}
transformed parameters {
  matrix[Ni, Nf] Lambda_mat = rep_matrix(0, Ni, Nf);

  {
    int pos = 0;
    for (i in 1:Ni) {
      for (j in 1:Nf) {
        if (Load_Pattern[i, j] == 1) {
          pos += 1;
          Lambda_mat[i, j] = lambda[pos];
        } else if (Load_Pattern[i, j] == -1) {
          pos += 1;
          Lambda_mat[i, j] = lambda[pos] * -1;
        }
      }
    }
  }
}
model {
  vector[NpNi] prec_s = precisions[i_ids];
  vector[NpNi] bb_a;
  vector[NpNi] bb_b;
  vector[NpNi - count_1 - count_max] lp_upper;
  vector[NpNi - count_1 - count_max] lp_lower;
  real lp_1;
  real lp_max;

  for (i in 1:NpNi) {
    bb_a[i] = inv_logit(means[i_ids[i]] + Lambda_mat[i_ids[i], ] * eta[p_ids[i]]) * prec_s[i];
  }
  bb_b = prec_s - bb_a;

  mu_m ~ student_t(3, 0, 1);
  sigma_m ~ student_t(3, 0, 1);
  means ~ normal(mu_m, sigma_m);

  lambda ~ lognormal(lambda_median, lambda_scale);
  eta ~ multi_normal_cholesky(rep_vector(0, Nf), R_chol);

  gamma_a ~ gamma(2, .1);
  gamma_b ~ gamma(2, .1);
  precisions ~ gamma(gamma_a, gamma_b);

  if (cuts[2] == 0) lp_1 = beta_binomial_lpmf(0 | N, bb_a[pos_1], bb_b[pos_1]);
  else lp_1 = beta_binomial_lcdf(cuts[2] | N, bb_a[pos_1], bb_b[pos_1]);

  lp_max = beta_binomial_lccdf(cuts[N_int] | N, bb_a[pos_max], bb_b[pos_max]);

  {
    int pos = 0;
    for (i in 1:NpNi) {
        if (n[i] > 1 && n[i] < N_int) {
        pos += 1;
        lp_upper[pos] = beta_binomial_lcdf(cuts[n[i] + 1] | N, bb_a[i], bb_b[i]);
        lp_lower[pos] = beta_binomial_lcdf(cuts[n[i]] | N, bb_a[i], bb_b[i]);
        }
    }
  }

  target += lp_1;
  target += lp_max;
  target += log_diff_exp(lp_upper, lp_lower);
}
generated quantities {
  vector[Ni] loc;
  vector[Ni] rho = rep_vector(0, Ni);
  vector[Ni] alpha_rel;
  vector[Ni] scale;
  real rho_m;
  real alpha_rel_m;
  matrix[Nf, Nf] R = multiply_lower_tri_self_transpose(R_chol);

  {
    vector[Ni] p = inv_logit(means);
    loc = N * p;
    rho = 1.0 ./ (precisions + 1);
    alpha_rel = rho ./ (rho + (1 - rho) / Ni);
    scale = sqrt(N * p .* (1 - p) .* (1 + (N - 1) * rho));
    rho_m = 1 / (gamma_a / gamma_b + 1);
    alpha_rel_m = rho_m / (rho_m + (1 - rho_m) / Ni);
  }
}
