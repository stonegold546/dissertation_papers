# REAL DATA ANALYSIS
# Implement a rounding strategy of people taking nearby numbers to round numbers

# Create own plot folder
plot_folder <- "~/Dropbox/dissertation/paper_3/plots/"

library(MASS)
library(rstan)
library(cmdstanr)
library(data.table)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

set_cmdstan_path("~/cmdstan/")

X <- read.csv("ICPSR_37106/DS0001/37106-0001-Data.tsv", sep = "\t")
colnames(X)
X <- X[, c(paste0("MH00", c(6, 7)), paste0("MH0", c(11, 12, 14, 16, 18, 20)))]
apply(X, 2, table)
X[X > 4] <- NA
colnames(X) <- paste0("x.", 1:8)

X.df.full <- reshape(X, paste0("x.", 1:8), direction = "long")
X.df.full <- X.df.full[complete.cases(X.df.full), ]
X.df.full$id.f <- as.integer(as.factor(X.df.full$id))
sum(table(X.df.full$id.f) == 8)
# [1] 6367

library(ggplot2)
library(ggforce)
library(scales)

X.df.dt <- as.data.table(X.df.full)
X.df.dt <- X.df.dt[, .N, list(time, x)]
X.df.dt[, time.r := abs(time - 9)]
X.df.dt[, x.f := factor(x, 1:4, c("< 1", "1 - 2", "3 - 4", "5 - 7"))]
X.df.dt[, x.f.r := factor(x, 4:1, rev(c("< 1", "1 - 2", "3 - 4", "5 - 7")))]
X.df.dt

# c(6, 7, 11, 12, 14, 16, 18, 20)
items <- paste0(1:8, ". ", c(
  "felt depressed", "everything an effort", "restless sleep",
  "was happy", "felt lonely", "enjoyed life", "felt sad", "not get 'going'"))

ggplot(X.df.dt, aes(time, N, col = x.f)) +
  geom_linerange(aes(ymin = 0, ymax = N), position = position_dodge(.6), size = 1) +
  theme_classic() +
  scale_y_continuous(trans = sqrt_trans(), breaks = c(0, 25, 200, 500, 1e3, seq(2e3, 6e3, 2e3))) +
  scale_x_continuous(breaks = 1:8, labels = items) +
  scale_color_manual(values = cbPalette) +
  theme(strip.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_line(size = .5),
        panel.grid.major.y = element_line(size = .5),
        panel.grid.minor.y = element_line(size = .25),
        axis.text.x = element_text(angle = 30, vjust = 1.1, hjust = 1),
        axis.ticks.x = element_blank(),
        legend.position = "top") +
  labs(x = "Items", y = "Frequency (square-root transformed axis)", shape = "Days", col = "Days")
ggsave(paste0(plot_folder, "real_freq_line.pdf"), width = 6.5, height = 4.5)

X.df.dt <- X.df.dt[order(time, x)]
X.df.dt[, N.t := ifelse(x == 4 & N <= 530, cumsum(N), cumsum(N) - N + N / 2), list(time)]
# X.df.dt[, N.t := cumsum(N) - N + N / 2, list(time)]
X.df.dt
X.df.dt[, N.p := N / sum(N), time]
X.df.dt[, items := items[time]]
X.df.dt

ggplot(X.df.dt, aes(reorder(items, time.r), N, fill = x.f.r)) +
  geom_bar(stat = "identity", width = .75, position = position_stack()) +
  geom_text(aes(y = N.t, label = paste0(number(N, 1, big.mark = ","))), # percent(N.p, 1), "\n", 
            col = 1, size = 2.5, angle = 15) +
  coord_flip() + theme_classic() +
  scale_y_continuous(breaks = c(seq(0, 5e3, 1e3), 6400), sec.axis = sec_axis(
    trans = ~ . / 6400, labels = percent_format(), breaks = seq(0, 1, .2),
    name = "Percentage of respondents")) +
  # scale_fill_manual(values = cbPalette[4:1]) +
  scale_fill_grey(start = .3, end = .9) +
  guides(fill = guide_legend(reverse = TRUE), col = guide_legend(reverse = TRUE)) +
  theme(legend.position = "top") +
  labs(x = "Items", y = "Number of respondents", fill = "Days in last week")
ggsave(paste0(plot_folder, "real_freq_stack.pdf"), width = 6.5, height = 4)

dat.list.mini <- list(
  N = 7, n_1d = X.df.dt$x, i_ids = X.df.dt$time, counts = X.df.dt$N)
dat.list.mini$Ni <- max(dat.list.mini$i_ids)
dat.list.mini$Nc <- length(dat.list.mini$n_1d)
dat.list.mini$N_int <- 4
(dat.list.mini$cuts <- matrix(c(0, 0, 1, 2, 3, 4, 5, 7), byrow = TRUE, ncol = 2))

int_bb <- cmdstan_model("stan_scripts/interval_bb_item.stan")
i.bb.fit <- int_bb$sample(
  data = dat.list.mini, seed = 12345, iter_warmup = 500,
  iter_sampling = 500, chains = 3, parallel_chains = 3)
i.bb.fit$cmdstan_diagnose()
i.bb.fit <- read_stan_csv(i.bb.fit$output_files())
print(i.bb.fit, c("log_lik", "scale", "p_lgt", "r_lgt"), digits_summary = 3, include = FALSE)

# Analyze subset of data first:

set.seed(12345)
X.df <- X.df.full[X.df.full$id %in% sample(unique(X.df.full$id), 1e3), ]
X.df$id.f <- as.integer(as.factor(X.df$id))

dat.list.int <- list(
  Np = max(X.df$id.f), Ni = max(X.df$time), N = 7, shape_r = 5,
  lambda_median = log(.5), lambda_scale = log(4 / .5) / qnorm(.99),
  neg_l = array(c(4, 6)))
dat.list.int$n_neg <- length(dat.list.int$neg_l)
dat.list.int$p_ids <- X.df$id.f
dat.list.int$i_ids <- X.df$time
dat.list.int$n_1d <- X.df$x
dat.list.int$NpNi <- length(dat.list.int$n_1d)
dat.list.int$N_int <- 4
# dat.list.int$cuts <- c(0, 0, 2, 4, 7)
(dat.list.int$cuts <- matrix(c(0, 0, 1, 2, 3, 4, 5, 7), byrow = TRUE, ncol = 2))
dat.list.int$count_1 <- sum(dat.list.int$n_1d == 1)
dat.list.int$count_max <- sum(dat.list.int$n_1d == 4)
dat.list.int$pos_1 <- which(dat.list.int$n_1d == 1)
dat.list.int$pos_max <- which(dat.list.int$n_1d == 4)

int_multi_bb <- cmdstan_model("stan_scripts/interval_multi_bb.stan")
i.mbb.fit <- int_multi_bb$sample(
  data = dat.list.int, seed = 12345, iter_warmup = 750, refresh = 50,
  iter_sampling = 750, chains = 3, parallel_chains = 3)
i.mbb.fit$cmdstan_diagnose()
i.mbb.fit <- read_stan_csv(i.mbb.fit$output_files())
print(i.mbb.fit, c("lambda_p", "eta"), digits_summary = 3, include = FALSE)
rowSums(get_elapsed_time(i.mbb.fit))
# chain:1 chain:2 chain:3 
# 277.383 316.166 301.697 
check_hmc_diagnostics(i.mbb.fit)
print(i.mbb.fit, c("sd_p_lgt", "m_p_lgt"), digits_summary = 4)
print(i.mbb.fit, c("p"), digits_summary = 3)
print(i.mbb.fit, c("rho"), digits_summary = 3)
print(i.mbb.fit, c("lambda"), digits_summary = 3)

int_multi_bb_lsm <- cmdstan_model("stan_scripts/interval_multi_bb_lsm.stan")
i.mbb.fit.lsm <- int_multi_bb_lsm$sample(
  data = dat.list.int, seed = 12345, iter_warmup = 750, refresh = 50,
  iter_sampling = 750, chains = 3, parallel_chains = 3)
i.mbb.fit.lsm$cmdstan_diagnose()
i.mbb.fit.lsm <- read_stan_csv(i.mbb.fit.lsm$output_files())
print(i.mbb.fit.lsm, c("Eta_pr"), digits_summary = 3, include = FALSE)
rowSums(get_elapsed_time(i.mbb.fit.lsm))
# chain:1 chain:2 chain:3 
# 426.861 420.597 451.490 
check_hmc_diagnostics(i.mbb.fit.lsm)
print(i.mbb.fit.lsm, c("r", "sd_p_lgt", "m_p_lgt", "sd_r_lgt", "m_r_lgt", "sigma_eta"),
      digits_summary = 4)
#              mean se_mean     sd    2.5%     25%     50%     75%   97.5% n_eff   Rhat
# r         -0.2195  0.0027 0.0897 -0.3898 -0.2795 -0.2205 -0.1590 -0.0423  1083 1.0005
# sd_p_lgt   2.3692  0.0106 0.5789  1.5200  1.9401  2.2939  2.7005  3.7123  2994 0.9994
# m_p_lgt   -1.7014  0.0135 0.7906 -3.1776 -2.2306 -1.7330 -1.1756 -0.0985  3442 0.9990
# sd_r_lgt   0.8652  0.0065 0.3022  0.4630  0.6485  0.7997  1.0082  1.6264  2142 0.9996
# m_r_lgt   -1.6046  0.0078 0.3447 -2.3050 -1.8074 -1.6164 -1.3959 -0.8700  1957 0.9996
# sigma_eta  1.2174  0.0042 0.1181  0.9917  1.1386  1.2173  1.2920  1.4626   786 1.0021
traceplot(i.mbb.fit.lsm, c("r", "sd_p_lgt", "m_p_lgt", "sd_r_lgt", "m_r_lgt", "sigma_eta"))
print(i.mbb.fit.lsm, c("p"), digits_summary = 3)
print(i.mbb.fit.lsm, c("rho"), digits_summary = 3)
print(i.mbb.fit.lsm, c("lambda"), digits_summary = 3)

orig.ids <- as.data.table(X.df)[, .N, list(id.f, id)][order(id.f), id]
table(outcome <- read.csv("ICPSR_37106/DS0001/37106-0001-Data.tsv", sep = "\t")[orig.ids, "SOCSOCPARTICIP"])
outcome[outcome > 1] <- NA
dat.list.int$outcome <- outcome
dat.list.int$outcome[is.na(dat.list.int$outcome)] <- -1
table(dat.list.int$outcome)

int_multi_bb_lsm_p1 <- cmdstan_model("stan_scripts/interval_multi_bb_lsm_1pred.stan")
i.mbb.fit.lsm.p1 <- int_multi_bb_lsm_p1$sample(
  data = dat.list.int, seed = 12345, iter_warmup = 750, refresh = 50,
  iter_sampling = 750, chains = 3, parallel_chains = 3)
i.mbb.fit.lsm.p1$cmdstan_diagnose()
i.mbb.fit.lsm.p1 <- read_stan_csv(i.mbb.fit.lsm.p1$output_files())
print(i.mbb.fit.lsm.p1, c("Eta_pr"), digits_summary = 3, include = FALSE)
rowSums(get_elapsed_time(i.mbb.fit.lsm.p1))
# chain:1 chain:2 chain:3 
# 543.953 526.171 562.719 
check_hmc_diagnostics(i.mbb.fit.lsm.p1)
print(i.mbb.fit.lsm.p1, digits_summary = 4,
      c("R[1,2]", "R[1,3]", "R[2,3]", "sd_p_lgt", "m_p_lgt", "sd_r_lgt", "m_r_lgt",
        "sigma_eta", "outcome_int", "prop", "coefs"))
#                mean se_mean     sd    2.5%     25%     50%     75%   97.5% n_eff   Rhat
# R[1,2]      -0.2193  0.0082 0.0846 -0.3778 -0.2765 -0.2229 -0.1648 -0.0414   107 1.0223
# R[1,3]      -0.1729  0.0015 0.0458 -0.2579 -0.2039 -0.1742 -0.1424 -0.0821   886 1.0032
# R[2,3]      -0.1468  0.0051 0.0767 -0.2958 -0.1989 -0.1473 -0.0940 -0.0033   228 1.0073
# sd_p_lgt     2.3826  0.0175 0.6281  1.5048  1.9320  2.2700  2.6996  3.9182  1288 0.9999
# m_p_lgt     -1.7310  0.0248 0.8139 -3.2648 -2.2428 -1.7807 -1.2409 -0.1027  1074 1.0005
# sd_r_lgt     0.8492  0.0095 0.2983  0.4581  0.6409  0.7890  0.9849  1.5961   981 0.9995
# m_r_lgt     -1.6191  0.0110 0.3288 -2.2781 -1.8202 -1.6240 -1.4200 -0.9694   888 1.0030
# sigma_eta    1.2191  0.0045 0.1182  0.9920  1.1405  1.2164  1.2951  1.4632   689 1.0030
# outcome_int -0.0136  0.0011 0.0383 -0.0874 -0.0398 -0.0138  0.0133  0.0587  1194 1.0002
# prop         0.4946  0.0004 0.0152  0.4652  0.4841  0.4945  0.5053  0.5234  1194 1.0002
# coefs[1]    -0.2163  0.0027 0.0536 -0.3191 -0.2520 -0.2177 -0.1799 -0.1131   406 1.0095
# coefs[2]    -0.1954  0.0053 0.0805 -0.3491 -0.2505 -0.1964 -0.1400 -0.0395   235 1.0054
traceplot(i.mbb.fit.lsm.p1, c("R[1,2]", "R[1,3]", "R[2,3]", "sd_p_lgt", "m_p_lgt",
                            "sd_r_lgt", "m_r_lgt", "sigma_eta", "outcome_int", "coefs"))
print(i.mbb.fit.lsm.p1, c("p"), digits_summary = 3)
print(i.mbb.fit.lsm.p1, c("rho"), digits_summary = 3)
print(i.mbb.fit.lsm.p1, c("lambda"), digits_summary = 3)

# Analyze full dataset:

dat.list.int.full <- list(
  Np = max(X.df.full$id.f), Ni = max(X.df.full$time), N = 7, shape_r = 2,
  lambda_median = log(.5), lambda_scale = log(4 / .5) / qnorm(.99),
  neg_l = array(c(4, 6)))
dat.list.int.full$n_neg <- length(dat.list.int.full$neg_l)
dat.list.int.full$p_ids <- X.df.full$id.f
dat.list.int.full$i_ids <- X.df.full$time
dat.list.int.full$n_1d <- X.df.full$x
dat.list.int.full$NpNi <- length(dat.list.int.full$n_1d)
dat.list.int.full$N_int <- 4
(dat.list.int.full$cuts <- matrix(c(0, 0, 1, 2, 3, 4, 5, 7), byrow = TRUE, ncol = 2))
dat.list.int.full$count_1 <- sum(dat.list.int.full$n_1d == 1)
dat.list.int.full$count_max <- sum(dat.list.int.full$n_1d == 4)
dat.list.int.full$pos_1 <- which(dat.list.int.full$n_1d == 1)
dat.list.int.full$pos_max <- which(dat.list.int.full$n_1d == 4)

int_multi_bb <- cmdstan_model("stan_scripts/interval_multi_bb.stan")
i.mbb.fit.full <- int_multi_bb$sample(
  data = dat.list.int.full, seed = 12345, iter_warmup = 750,
  iter_sampling = 750, chains = 3, parallel_chains = 3, refresh = 50)
i.mbb.fit.full$cmdstan_diagnose()
i.mbb.fit.full <- read_stan_csv(i.mbb.fit.full$output_files())
print(i.mbb.fit.full, c("eta"), digits_summary = 3, include = FALSE)
rowSums(get_elapsed_time(i.mbb.fit.full))
#  chain:1  chain:2  chain:3 
# 3030.693 2965.883 3012.150 
check_hmc_diagnostics(i.mbb.fit.full)
print(i.mbb.fit.full, c("sd_p_lgt", "m_p_lgt"), digits_summary = 4)
print(i.mbb.fit.full, c("p"), digits_summary = 3)
print(i.mbb.fit.full, c("rho"), digits_summary = 3)
print(i.mbb.fit.full, c("lambda"), digits_summary = 3)

int_multi_bb_lsm <- cmdstan_model("stan_scripts/interval_multi_bb_lsm.stan")
i.mbb.fit.lsm.full <- int_multi_bb_lsm$sample(
  data = dat.list.int.full, seed = 12345, iter_warmup = 750,
  iter_sampling = 750, chains = 3, parallel_chains = 3, refresh = 50)
i.mbb.fit.lsm.full$cmdstan_diagnose()
i.mbb.fit.lsm.full <- read_stan_csv(i.mbb.fit.lsm.full$output_files())
print(i.mbb.fit.lsm.full, c("Eta_pr"), digits_summary = 3, include = FALSE)
rowSums(get_elapsed_time(i.mbb.fit.lsm.full))
# chain:1 chain:2 chain:3 
# 3926.81 3951.64 3820.41 
check_hmc_diagnostics(i.mbb.fit.lsm.full)
print(i.mbb.fit.lsm.full, c("r", "sd_p_lgt", "m_p_lgt", "sd_r_lgt", "m_r_lgt", "sigma_eta"),
      digits_summary = 4)
#              mean se_mean     sd    2.5%     25%     50%     75%   97.5% n_eff   Rhat
# r         -0.1624  0.0015 0.0371 -0.2347 -0.1875 -0.1623 -0.1366 -0.0908   623 1.0038
# sd_p_lgt   2.3955  0.0120 0.6038  1.4900  1.9592  2.3002  2.7262  3.8669  2541 0.9992
# m_p_lgt   -1.7427  0.0153 0.7751 -3.2436 -2.2443 -1.7684 -1.2446 -0.1469  2576 0.9999
# sd_r_lgt   0.8048  0.0072 0.2621  0.4749  0.6270  0.7486  0.9181  1.4933  1343 1.0001
# m_r_lgt   -1.7377  0.0076 0.3004 -2.3363 -1.9208 -1.7473 -1.5680 -1.1089  1579 1.0001
# sigma_eta  1.3445  0.0021 0.0504  1.2524  1.3089  1.3441  1.3799  1.4457   591 1.0064
traceplot(i.mbb.fit.lsm.full, c("r", "sd_p_lgt", "m_p_lgt", "sd_r_lgt", "m_r_lgt", "sigma_eta"))
print(i.mbb.fit.lsm.full, c("p"), digits_summary = 3)
print(i.mbb.fit.lsm.full, c("rho"), digits_summary = 3)
print(i.mbb.fit.lsm.full, c("lambda"), digits_summary = 3)

round(plogis(qnorm(c(.025, .25, .5, .75, .975), -1.7377, 1.3445)), 3)
# [1] 0.012 0.066 0.150 0.303 0.710

orig.ids <- as.data.table(X.df.full)[, .N, list(id.f, id)][order(id.f), id]
table(outcome <- read.csv("ICPSR_37106/DS0001/37106-0001-Data.tsv", sep = "\t")[orig.ids, "SOCSOCPARTICIP"])
outcome[outcome > 1] <- NA
dat.list.int.full$outcome <- outcome
dat.list.int.full$outcome[is.na(dat.list.int.full$outcome)] <- -1
table(dat.list.int.full$outcome)

int_multi_bb_lsm_p1 <- cmdstan_model("stan_scripts/interval_multi_bb_lsm_1pred.stan")
i.mbb.fit.lsm.p1.full <- int_multi_bb_lsm_p1$sample(
  data = dat.list.int.full, seed = 12345, iter_warmup = 750, refresh = 50,
  iter_sampling = 750, chains = 3, parallel_chains = 3)
i.mbb.fit.lsm.p1.full$cmdstan_diagnose()
i.mbb.fit.lsm.p1.full <- read_stan_csv(i.mbb.fit.lsm.p1.full$output_files())
print(i.mbb.fit.lsm.p1.full, c("Eta_pr"), digits_summary = 3, include = FALSE)
rowSums(get_elapsed_time(i.mbb.fit.lsm.p1.full))
# chain:1 chain:2 chain:3 
# 6693.48 6714.76 6541.61 
check_hmc_diagnostics(i.mbb.fit.lsm.p1.full)
print(i.mbb.fit.lsm.p1.full, digits_summary = 4,
      c("R[1,2]", "R[1,3]", "R[2,3]", "sd_p_lgt", "m_p_lgt", "sd_r_lgt", "m_r_lgt",
        "sigma_eta", "outcome_int", "prop", "coefs"))
#                mean se_mean     sd    2.5%     25%     50%     75%   97.5% n_eff   Rhat
# R[1,2]      -0.1625  0.0024 0.0379 -0.2364 -0.1883 -0.1630 -0.1367 -0.0871   248 1.0066
# R[1,3]      -0.1494  0.0005 0.0181 -0.1848 -0.1615 -0.1491 -0.1377 -0.1127  1467 1.0004
# R[2,3]      -0.0936  0.0013 0.0283 -0.1497 -0.1130 -0.0936 -0.0744 -0.0389   486 1.0073
# sd_p_lgt     2.4214  0.0105 0.5898  1.5158  1.9933  2.3312  2.7690  3.7900  3150 0.9989
# m_p_lgt     -1.7540  0.0153 0.8377 -3.3414 -2.2951 -1.7734 -1.2574  0.0274  2985 1.0000
# sd_r_lgt     0.8020  0.0056 0.2652  0.4486  0.6152  0.7472  0.9228  1.4904  2255 0.9996
# m_r_lgt     -1.7361  0.0069 0.2986 -2.3119 -1.9170 -1.7389 -1.5567 -1.1144  1897 0.9989
# sigma_eta    1.3449  0.0019 0.0506  1.2482  1.3100  1.3457  1.3794  1.4455   710 1.0053
# outcome_int  0.0058  0.0003 0.0163 -0.0248 -0.0055  0.0053  0.0165  0.0375  3370 0.9993
# prop         0.5023  0.0001 0.0065  0.4901  0.4978  0.5021  0.5066  0.5149  3370 0.9993
# coefs[1]    -0.1692  0.0007 0.0197 -0.2083 -0.1822 -0.1687 -0.1560 -0.1315   901 1.0019
# coefs[2]    -0.1213  0.0013 0.0290 -0.1783 -0.1409 -0.1210 -0.1023 -0.0647   472 1.0099
traceplot(i.mbb.fit.lsm.p1.full, c("R[1,2]", "R[1,3]", "R[2,3]", "sd_p_lgt", "m_p_lgt",
                                   "sd_r_lgt", "m_r_lgt", "sigma_eta", "outcome_int", "coefs"))
print(i.mbb.fit.lsm.p1.full, c("p"), digits_summary = 3)
print(i.mbb.fit.lsm.p1.full, c("rho"), digits_summary = 3)
print(i.mbb.fit.lsm.p1.full, c("lambda"), digits_summary = 3)

round(apply(c.s <- (t(apply(as.data.frame(i.mbb.fit.lsm.p1.full, "R"), 1, function (rs) {
  R <- matrix(as.numeric(rs), 3)
  R.inv <- solve(R[1:2, 1:2])
  r.yx <- R[1:2, 3]
  beta.s <- R.inv %*% r.yx
  c(coefs = beta.s, r2 = r.yx %*% beta.s)
}))), 2, quantile, probs = c(.025, .5, .975)), 4)
#        coefs1  coefs2     r2
# 2.5%  -0.2083 -0.1783 0.0218
# 50%   -0.1687 -0.1210 0.0370
# 97.5% -0.1315 -0.0647 0.0573
round(colMeans(c.s), 4)
#  coefs1  coefs2      r2 
# -0.1692 -0.1213  0.0378 

summary(lavaan::sem(
  "F =~ x.1 + x.2 + x.3 + x.4 + x.5 + x.6 + x.7 + x.8
  outcome ~ F",
  data = cbind(X[orig.ids, ], outcome = outcome),
  std.lv = TRUE, ordered = c(paste0("x.", 1:8), "outcome")),
  rsquare = TRUE, standardize = TRUE, fit.measures = TRUE)

dim(theta.s <- t(as.data.frame(i.mbb.fit.lsm.p1.full, "Eta_pr")[, seq(1, 12790, 2)]))
dim(eta.s <- t(as.data.frame(i.mbb.fit.lsm.p1.full, "Eta_pr")[, seq(2, 12790, 2)]))
int.s <- as.data.frame(i.mbb.fit.lsm.p1.full, "outcome_int")
c.s <- as.data.frame(i.mbb.fit.lsm.p1.full, "coefs")

hist(plogis(-colMeans(int.s) + colMeans(c.s)[1] * rowMeans(theta.s) + colMeans(c.s)[2] * rowMeans(eta.s)))
plot(plogis(-colMeans(int.s) + colMeans(c.s)[1] * rowMeans(theta.s) + colMeans(c.s)[2] * rowMeans(eta.s))[
  dat.list.int.full$outcome >= 0],
  dat.list.int.full$outcome[dat.list.int.full$outcome >= 0])
print(psych::describeBy(plogis(
  -colMeans(int.s) + colMeans(c.s)[1] * rowMeans(theta.s) + colMeans(c.s)[2] * rowMeans(eta.s))[
    dat.list.int.full$outcome >= 0],
  dat.list.int.full$outcome[dat.list.int.full$outcome >= 0]), digits = 5)

round(t(apply(c.dist <- t(sapply(1:ncol(theta.s), function (i) {
  z.s <- dnorm(- int.s[i, 1] + c.s[i, 1] * theta.s[, i] + c.s[i, 2] * eta.s[, i])
  c(mean(z.s * c.s[i, 1]), mean(z.s * c.s[i, 2]))
})), 2, quantile, probs = c(.025, .5, .975))), 5)
#          2.5%      50%    97.5%
# [1,] -0.08025 -0.06586 -0.05169
# [2,] -0.06888 -0.04722 -0.02547

rm(theta.s, eta.s, int.s, c.s)
gc()

(c.dist.dt <- as.data.table(c.dist))
setnames(c.dist.dt, c("V1", "V2"), c("location", "ICC"))
c.dist.dt

c.dist.dt <- as.data.table(apply(c.dist.dt, 2, quantile, seq(0, 1, .001)))
c.dist.dt
c.dist.dt$probs <- seq(0, 1, .001)
c.dist.dt
round(c.dist.dt[probs %in% c(.025, .5, .975)], 3)
#    location    ICC probs
# 1:   -0.080 -0.069 0.025
# 2:   -0.066 -0.047 0.500
# 3:   -0.052 -0.025 0.975

# (ggplot(c.dist.dt, aes(location, probs)) + geom_line() +
#     scale_x_continuous(breaks = c.dist.dt[round(probs, 4) %in% c(.025, .25, .5, .75, .975), location],
#                        labels = number_format(.001)) +
#     scale_y_continuous(breaks = c.dist.dt[round(probs, 4) %in% c(.025, .25, .5, .75, .975), probs],
#                        labels = percent_format(.1)) +
#     theme_bw() + theme(panel.grid.minor.x = element_blank()) +
#     labs(y = "Prob(coefficient > threshold)", x = "location coefficient")) /
#   (ggplot(c.dist.dt, aes(ICC, probs)) + geom_line() +
#      scale_x_continuous(breaks = c.dist.dt[round(probs, 4) %in% c(.025, .25, .5, .75, .975), ICC],
#                         labels = number_format(.001)) +
#      scale_y_continuous(breaks = c.dist.dt[round(probs, 4) %in% c(.025, .25, .5, .75, .975), probs],
#                         labels = percent_format(.1)) +
#      theme_bw() + theme(panel.grid.minor.x = element_blank()) +
#      labs(y = "Prob(coefficient > threshold)", x = "ICC coefficient"))

i.df <- as.data.frame(summary(
  i.bb.fit, c("p", "rho"), probs = c(.025, .5, .975))$summary)
i.df$params <- rownames(i.df)

(mu.s <- i.df$mean[1:8])
(prec.s <- exp(-qlogis(i.df$mean[9:16])))
(alpha.s <- mu.s * prec.s)
(beta.s <- prec.s - alpha.s)
(i.df.qty <- as.data.table(t(sapply(1:8, function (i) TailRank::qbb(
  c(.025, .25, .5, .75, .975), 7, alpha.s[i], beta.s[i])))))
colnames(i.df.qty) <- c("lll", "ll", "md", "ul", "ull")
i.df.qty
i.df.qty[, items := items]
i.df.qty[, time.r := 8:1]
i.df.qty$mean <- mu.s * 7
i.df.qty

ggplot(i.df.qty, aes(reorder(items, time.r), md)) +
  coord_flip() + geom_point(size = 4, shape = 4) +
  geom_text(aes(y = mean, label = paste0(number(mean, .1), " days")), vjust = -.65,
            data = i.df.qty[1, ]) +
  geom_text(aes(y = mean, label = number(mean, .1)), vjust = -.65,
            data = i.df.qty[-1, ]) +
  theme_bw() +
  geom_point(aes(y = mean), shape = 1, size = 3) +
  geom_linerange(aes(ymin = ll, ymax = ul), size = 4, alpha = .25) +
  geom_linerange(aes(ymin = lll, ymax = ull)) +
  scale_y_continuous(breaks = 0:7) +
  labs(y = "Average with middle 50% and 95%",
       subtitle = "Median = X, mean = O number of days") +
  theme(panel.border = element_blank(), axis.ticks = element_blank(),
        axis.title.y = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank())
ggsave(paste0(plot_folder, "base_descriptives.pdf"), width = 6.5, height = 3.5)

i.fit.df <- as.data.frame(summary(
  i.mbb.fit, c("lambda", "p", "rho"), probs = c(.025, .5, .975))$summary)
i.fit.df$params <- rownames(i.fit.df)
i.lsm.fit.df <- as.data.frame(summary(
  i.mbb.fit.lsm, c("r", "sigma_eta", "lambda", "p", "rho"), probs = c(.025, .5, .975))$summary)
i.lsm.fit.df$params <- rownames(i.lsm.fit.df)
i.full.fit.df <- as.data.frame(summary(
  i.mbb.fit.full, c("lambda", "p", "rho"), probs = c(.025, .5, .975))$summary)
i.full.fit.df$params <- rownames(i.full.fit.df)
i.lsm.full.fit.df <- as.data.frame(summary(
  i.mbb.fit.lsm.full, c("r", "sigma_eta", "lambda", "p", "rho"), probs = c(.025, .5, .975))$summary)
i.lsm.full.fit.df$params <- rownames(i.lsm.full.fit.df)
i.lsm.full.fit.p1.df <- as.data.frame(summary(
  i.mbb.fit.lsm.p1.full, c("R[1,2]", "R[1,3]", "R[2,3]", "coefs", "sigma_eta", "lambda", "p", "rho"),
  probs = c(.025, .5, .975))$summary)
i.lsm.full.fit.p1.df$params <- rownames(i.lsm.full.fit.p1.df)

comb.df <- rbindlist(list(
  i.df, i.fit.df, i.lsm.fit.df, i.full.fit.df, i.lsm.full.fit.df, i.lsm.full.fit.p1.df), idcol = "fit")
comb.df

comb.df[, class := gsub("\\[\\d+\\]", "", params)]
comb.df[, class := ifelse(class == "lambda", "loading", class)]
comb.df[, class := ifelse(class == "p", "mean X 7 days", class)]
comb.df[, class := ifelse(class == "rho", "ICC", class)]
comb.df[, class := ifelse(class == "sigma_eta", "load_icc", class)]
comb.df[, item := as.integer(gsub("[a-z]+|\\[|\\]|\\_", "", params))]
setnames(comb.df, c("50%", "2.5%", "97.5%"), c("median", "ll", "ul"))
comb.df

comb.df[, pos := ifelse(
  class == "loading", as.integer(item), ifelse(
    class == "mean", as.integer(item) + 8, ifelse(
      class == "ICC", as.integer(item) + 16, 17)))]

comb.df[, fit.t := factor(
  fit, labels = c("Baseline", "Location (n = 1000)", "Location-ICC (n = 1000)",
                  "Location", "Location-ICC", "Full"))]
comb.df[, neg := item %in% c(4, 6)]
comb.df[class == "mean X 7 days", median := median * 7]
comb.df[class == "mean X 7 days", ll := ll * 7]
comb.df[class == "mean X 7 days", ul := ul * 7]
# comb.df[class == "mean" & neg, median := 1 - median]
# comb.df[class == "mean" & neg, ll := 1 - ll]
# comb.df[class == "mean" & neg, ul := 1 - ul]

library(ggplot2)
library(ggforce)
library(scales)
library(patchwork)

items <- paste0(1:8, ". ", c(
  "felt depressed", "everything an effort", "restless sleep",
  "was happy", "felt lonely", "enjoyed life", "felt sad", "not get 'going'"))

comb.df[, .N, list(fit, fit.t)]
ggplot(comb.df[class %in% c("ICC", "loading", "mean X 7 days") & fit %in% c(1, 4:6)],
       aes(item, abs(median))) +
  facet_col(~ reorder(class, pos), scales = "free_y", strip.position = "right") +
  geom_rect(xmin = 3.5, xmax = 4.5, ymin = -5, ymax = 50, fill = "white",
            size = .125, alpha = .025, col = "#999999") +
  geom_rect(xmin = 5.5, xmax = 6.5, ymin = -5, ymax = 50, fill = "white",
            size = .125, alpha = .025, col = "#999999") +
  geom_point(aes(shape = fit.t), position = position_dodge(.5)) +
  geom_linerange(aes(ymin = abs(ll), ymax = abs(ul), fill = fit.t),
                 position = position_dodge(.5), size = .3) +
  # scale_color_manual(values = cbPalette) +
  scale_shape_manual(values = c(1, 5, 4, 8)) +
  scale_x_continuous(breaks = 1:8, labels = items) + theme_bw() +
  theme(axis.text.x = element_text(angle = 30, vjust = 1.1, hjust = 1),
        legend.position = "top", axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank()) +
  labs(x = "Item", y = "Median with 95% quantile intervals", fill = "", col = "Model", shape = "Model")
ggsave(paste0(plot_folder, "real_params.pdf"), width = 6.5, height = 5)

Eta.pr.s <- as.data.frame(i.mbb.fit.lsm.full, "Eta_pr")
Rchol <- as.data.frame(i.mbb.fit.lsm.full, "R_chol")
Eta.s <- list()
for (i in 1:nrow(Rchol)) {
  Eta.s[[i]] <- t(matrix(as.numeric(Rchol[i, ]), 2) %*%
                    matrix(as.numeric(Eta.pr.s[i, ]), 2))
}
rm(Eta.pr.s)
rm(Rchol)

p.s <- as.data.frame(i.mbb.fit.lsm.full, "p_lgt")
lambda.s <- as.data.frame(i.mbb.fit.lsm.full, "lambda")
prob.s <- sapply(1:nrow(p.s), function (i) {
  plogis(
    as.numeric(p.s[i, ])[dat.list.int.full$i_ids] +
      as.numeric(lambda.s[i, ])[dat.list.int.full$i_ids] *
      Eta.s[[i]][dat.list.int.full$p_ids, 1])
})
dim(prob.s)
rm(p.s, lambda.s)

r.s <- as.data.frame(i.mbb.fit.lsm.full, "r_lgt")
load_icc <- as.data.frame(i.mbb.fit.lsm.full, "sigma_eta")
prec.s <- sapply(1:nrow(r.s), function (i) {
  exp(
    -as.numeric(r.s[i, ])[dat.list.int.full$i_ids] -
      as.numeric(load_icc[i, ]) * Eta.s[[i]][dat.list.int.full$p_ids, 2])
})
dim(prec.s)
rm(r.s, load_icc)

set.seed(12345)
N.tab <- matrix(nrow = ncol(prob.s) * dat.list.int.full$Ni, ncol = 9)
pos <- 0
for (j in 1:ncol(prob.s)) {
  for (i in 1:dat.list.int.full$Ni) {
    idxs <- dat.list.int.full$i_ids == i
    len <- sum(idxs)
    p <- prob.s[idxs, j]
    prec <- prec.s[idxs, j]
    p <- rbeta(len, p * prec, (1 - p) * prec)
    counts <- rbinom(len, 7, p)
    pos <- pos + 1
    N.tab[pos, ] <- c(i, sapply(0:7, function (i) sum(counts == i)))
  }
}
rm(j, p, counts, pos, idxs, len)
rm(prob.s, prec.s)
# rm(Eta.s)
gc()

N.tab <- as.data.table(N.tab)
colnames(N.tab) <- c("item", paste0("count.", 0:7))
N.tab

N.tab.l <- reshape(N.tab, paste0("count.", 0:7), direction = "long")

org.dat <- as.data.table(X.df.full)[, .N, list(item = time, x)]
org.dat <- as.data.table(merge(data.frame(x = c(1, 2, 2, 3, 3, 4, 4, 4), time = 0:7), org.dat))
org.dat <- org.dat[order(item)]
org.dat[, N := N / ifelse(time == 0, 1, ifelse(time <= 4, 2, 3))]
org.dat

ggplot(N.tab.l, aes(time, count)) +
  geom_segment(aes(x = time - .35, xend = time + .35, y = count, yend = count),
               alpha = .5, size = .25, col = cbPalette[1]) +
  # scale_x_reverse() +
  scale_y_continuous(trans = sqrt_trans(), breaks = c(20, 200, 500, 1000, seq(2e3, 8e3, 2e3))) +
  facet_wrap(~ reorder(paste0("item ", item), item), ncol = 4, strip.position = "top") + theme_bw() +
  # coord_flip() +
  geom_segment(aes(x = time - .5, xend = time + .5, y = N, yend = N), data = org.dat) +
  theme(strip.background = element_blank(),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  labs(x = "Days", y = "Frequency (square root transformed axis)")
ggsave(paste0(plot_folder, "real_posterior_samples.pdf"), width = 6.5, height = 4)

orig.ids <- as.data.table(X.df.full)[, .N, list(id.f, id)][order(id.f), id]

lt.mean <- as.data.table(lapply(Eta.s, function (x) x[, 1]))
lt.mean <- as.data.table(t(apply(lt.mean, 1, quantile, probs = c(.5, .25, .75))))
setnames(lt.mean, c("50%", "25%", "75%"), c("median_m", "ll_m", "ul_m"))
lt.mean$score <- as.matrix(X[orig.ids, ]) %*% c(rep(1, 3), -1, 1, -1, rep(1, 2))
lt.mean$mdn_reduced <- apply(as.data.frame(i.mbb.fit.full, "eta"), 2, median)

plot(mdn_reduced ~ median_m, lt.mean)
cor(lt.mean$mdn_reduced, lt.mean$median_m)
# IMPORTANT: Very high correlation in latents
cor(lt.mean$mdn_reduced, lt.mean$median_m, method = "sp")
cor(lt.mean, use = "p")
cor(lt.mean, use = "p", method = "sp")

ggplot(lt.mean, aes(score, median_m)) +
  # geom_ribbon(aes(ymin = ll_m, ymax = ul_m), alpha = .5, size = .5, col = cbPalette[1]) +
  geom_linerange(aes(ymin = ll_m, ymax = ul_m), alpha = .5, size = .5, col = cbPalette[1]) +
  geom_jitter(height = 0, width = .1, shape = 4) + theme_bw()

lt.icc <- as.data.table(lapply(Eta.s, function (x) x[, 2]))
lt.icc <- as.data.table(t(apply(lt.icc, 1, quantile, probs = c(.5, .25, .75))))
setnames(lt.icc, c("50%", "25%", "75%"), c("median_icc", "ll_icc", "ul_icc"))
lt.mean <- cbind(lt.mean, lt.icc)
rm(lt.icc)

ggplot(lt.mean, aes(score, median_icc)) +
  geom_linerange(aes(ymin = ll_icc, ymax = ul_icc), alpha = .5, size = .5, col = cbPalette[1]) +
  geom_jitter(height = 0, width = .1, shape = 4) + theme_bw()

ggplot(lt.mean, aes(median_m, median_icc)) +
  geom_point(shape = 4) + theme_bw()

