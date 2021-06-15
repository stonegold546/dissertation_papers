data {
  int Ni;
  int Nc;
  int n_1d[Nc];
  int i_ids[Nc];
  int N;
  int N_int;
  int cuts[N_int, 2];
  int counts[Nc];
}
parameters {
  vector[Ni] p_lgt;
  vector[Ni] r_lgt;
}
model {
  p_lgt ~ student_t(3, 0, 2.5);
  r_lgt ~ student_t(3, 0, 2.5);

  {
    vector[Nc] p_s = inv_logit(p_lgt)[i_ids];
    vector[Nc] disp_s = exp(-r_lgt)[i_ids];
    vector[Nc] bb_a = p_s .* disp_s;
    vector[Nc] bb_b = disp_s - bb_a;

    {
      real lp_upper;
      real lp_lower;
      for (i in 1:Nc) {
        if (n_1d[i] == 1) {
          if (cuts[1, 2] == 0) target += beta_binomial_lpmf(0 | N, bb_a[i], bb_b[i]) * counts[i];
          else target += beta_binomial_lcdf(cuts[1, 2] | N, bb_a[i], bb_b[i]) * counts[i];
        } else if (n_1d[i] == N_int) {
          target += beta_binomial_lccdf(cuts[N_int, 1] - 1 | N, bb_a[i], bb_b[i]) * counts[i];
        } else {
          lp_upper = beta_binomial_lcdf(cuts[n_1d[i], 2] | N, bb_a[i], bb_b[i]);
          lp_lower = beta_binomial_lcdf(cuts[n_1d[i], 1] - 1 | N, bb_a[i], bb_b[i]);
          target += log_diff_exp(lp_upper, lp_lower) * counts[i];
        }
      }
    }
  }
}
generated quantities {
  vector[Ni] p = inv_logit(p_lgt);
  vector[Ni] rho = inv_logit(r_lgt);
  vector[Ni] scale = sqrt(N * p .* (1 - p) .* (1 + (N - 1) * rho));
  vector[Nc] log_lik;

  {
    vector[Nc] p_s = inv_logit(p_lgt)[i_ids];
    vector[Nc] disp_s = exp(-r_lgt)[i_ids];
    vector[Nc] bb_a = p_s .* disp_s;
    vector[Nc] bb_b = disp_s - bb_a;

    {
      real lp_upper;
      real lp_lower;
      for (i in 1:Nc) {
        if (n_1d[i] == 1) {
          if (cuts[1, 2] == 0) log_lik[i] = beta_binomial_lpmf(0 | N, bb_a[i], bb_b[i]) * counts[i];
          else log_lik[i] = beta_binomial_lcdf(cuts[1, 2] | N, bb_a[i], bb_b[i]) * counts[i];
        } else if (n_1d[i] == N_int) {
          log_lik[i] = beta_binomial_lccdf(cuts[N_int, 1] - 1 | N, bb_a[i], bb_b[i]) * counts[i];
        } else {
          lp_upper = beta_binomial_lcdf(cuts[n_1d[i], 2] | N, bb_a[i], bb_b[i]);
          lp_lower = beta_binomial_lcdf(cuts[n_1d[i], 1] - 1 | N, bb_a[i], bb_b[i]);
          log_lik[i] = log_diff_exp(lp_upper, lp_lower) * counts[i];
        }
      }
    }
  }
}
