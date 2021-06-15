# SIMULATION DATA ANALYSIS # 2

plot_folder <- "~/Dropbox/dissertation/paper_3/plots/"

library(MASS)
library(rstan)
library(cmdstanr)
library(data.table)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

set_cmdstan_path("~/cmdstan/")

beta_bin_lsm <- cmdstan_model("stan_scripts/multi_beta_bin_lsm_multidim.stan")
beta_bin_lsm_interval <- cmdstan_model("stan_scripts/interval_multi_bb_lsm_multidim.stan")

set.seed(12345)
params <- data.frame(
  means = sample(c(.05, .05, .1, .1, .2, .9), 6),
  iccs = sample(c(.2, .2, .2, .3, .3, .5), 6),
  lambdas = 1)
params
params$precs <- exp(-qlogis(params$iccs))
params$lambdas[params$means == .9] <- params$lambdas[params$means == .9] * -1
params

n.i <- nrow(params)

(Lambda <- matrix(c(params$lambdas[1:3], rep(0, 6), params$lambdas[4:6]), ncol = 2))

(R <- matrix(.5, 3, 3) + .5 * diag(3))
R[3, 1] <- R[1, 3] <- .35
R[3, 2] <- R[2, 3] <- .25
R

set.seed(12345)
cov(Eta <- mvrnorm(2e3, rep(0, 3), R, empirical = TRUE))
plot(Eta)
cov(M1 <- t(apply(Eta, 1, function (x) { Lambda %*% x[1:2] + qlogis(params$means) })))
plogis(colMeans(M1))
cov(M2 <- t(matrix(rep(1, n.i)) %*% Eta[, 3] + log(params$precs)))
exp(colMeans(M2))

cov(cbind(M1, M2))

set.seed(12345)
N <- matrix(nrow = nrow(M1), ncol = ncol(M1))
for (j in 1:ncol(M1)) {
  p <- plogis(M1[, j])
  prec <- exp(M2[, j])
  p <- rbeta(nrow(N), p * prec, (1 - p) * prec)
  N[, j] <- rbinom(nrow(N), 7, p)
}
rm(j, p)

plot(as.data.frame(N))

cor(as.data.frame(N))

apply(N, 2, median)
apply(N, 2, mean)
apply(N, 2, sd)

N.df <- as.data.frame(N)
colnames(N.df) <- paste0("n.", 1:ncol(N.df))
N.df <- na.omit(reshape(N.df, direction = "long", 1:ncol(N)))

N.df.dt <- as.data.table(N.df)
(N.df.dt <- N.df.dt[, .N, list(time, n)])
params$time <- 1:nrow(params)
N.df.dt <- merge(N.df.dt, params)
N.df.dt[, means.t := paste0("mean = ", means)]
N.df.dt[, iccs.t := paste0("ICC = ", iccs)]
N.df.dt[, factor := factor(ifelse(time < 4, 1, 2))]
N.df.dt

library(scales)

ggplot(N.df.dt, aes(n, N)) +
  geom_bar(col = 1, stat = "identity", position = position_dodge(.5)) +
  facet_wrap(~ reorder(paste0("item ", time, "\n", means.t, "\n", iccs.t), time), nrow = 1) +
  theme_bw() + theme(legend.position = "top") +
  labs(x = "Counts", y = "Frequency") +
  scale_fill_manual(values = cbPalette)
ggsave(paste0(plot_folder, "sim_item_dist_2.pdf"), width = 6.5, height = 3.5)

dat.list <- list(
  Np = nrow(N), Ni = ncol(N), NpNi = nrow(N.df), N = 7,
  n_1d = N.df$n, p_ids = N.df$id, i_ids = N.df$time, shape_r = 2,
  lambda_median = log(.5), lambda_scale = log(4 / .5) / qnorm(.99),
  Nf = ncol(Lambda), Nl = sum(Lambda != 0),
  Load_Pattern = (Lambda > 0) + -1 * (Lambda < 0))
# saveRDS(dat.list, "sim_dat_res/sim_dat_list_md.rds")
# dat.list <- readRDS("sim_dat_res/sim_dat_list_md.rds")

multi_beta_bin_lsm <- cmdstan_model("stan_scripts/multi_beta_bin_lsm_multidim.stan")
mbb.lsm.fit <- multi_beta_bin_lsm$sample(
  data = dat.list, seed = 12345, iter_warmup = 750,
  iter_sampling = 750, chains = 3, parallel_chains = 3)
mbb.lsm.fit$cmdstan_diagnose()
mbb.lsm.fit <- read_stan_csv(mbb.lsm.fit$output_files())
print(mbb.lsm.fit, c("Eta", "R_chol"), digits_summary = 3, include = FALSE)
# saveRDS(mbb.lsm.fit, "sim_dat_res/mbb.fit.lsm.md.rds")
# mbb.lsm.fit <- readRDS("sim_dat_res/mbb.fit.lsm.md.rds")
rowSums(get_elapsed_time(mbb.lsm.fit))
# chain:1 chain:2 chain:3 
# 412.848 390.260 400.484 

params$means
print(mbb.lsm.fit, c("p"), digits_summary = 3)

params$iccs
print(mbb.lsm.fit, c("rho"), digits_summary = 3)

L.idxs <- paste0(
  "Lambda_mat[", apply(which(Lambda != 0, arr.ind = TRUE), 1, paste0, collapse = ","), "]")
params$lambda
print(mbb.lsm.fit, L.idxs, digits_summary = 3)

R
print(mbb.lsm.fit, c("R[1,2]", "R[1,3]", "R[2,3]", "sigma_eta"), digits_summary = 3,
      probs = c(.025, .05, .25, .5, .75, .95, .975))
#             mean se_mean    sd   2.5%     5%    25%    50%    75%    95%  97.5% n_eff  Rhat
# R[1,2]     0.521   0.003 0.045  0.429  0.445  0.493  0.522  0.553  0.591  0.603   171 1.011
# R[1,3]    -0.291   0.006 0.069 -0.428 -0.406 -0.335 -0.291 -0.244 -0.177 -0.155   132 1.015
# R[2,3]    -0.244   0.005 0.068 -0.373 -0.352 -0.291 -0.248 -0.198 -0.128 -0.105   177 1.027
# sigma_eta  1.062   0.003 0.066  0.938  0.952  1.018  1.062  1.106  1.172  1.191   649 1.007

table(N.df$q <- cut(dat.list$n_1d, c(-1, 0, 2, 4, 7)))
table(N.df$q.int <- as.integer(N.df$q))

dat.list.int <- dat.list
dat.list.int$n_1d <- N.df$q.int
dat.list.int$N_int <- 4
# dat.list.int$cuts <- c(0, 0, 2, 4, 7)
(dat.list.int$cuts <- matrix(c(0, 0, 1, 2, 3, 4, 5, 7), byrow = TRUE, ncol = 2))
dat.list.int$count_1 <- sum(dat.list.int$n_1d == 1)
dat.list.int$count_max <- sum(dat.list.int$n_1d == 4)
dat.list.int$pos_1 <- which(dat.list.int$n_1d == 1)
dat.list.int$pos_max <- which(dat.list.int$n_1d == 4)
# dat.list.int$shape_r <- 3
# dat.list.int$lambda_scale <- log(3 / .5) / qnorm(.99)

dat.list.int$cuts
table(dat.list$n_1d, dat.list.int$n_1d)

# 101
int_multi_bb_lsm <- cmdstan_model("stan_scripts/interval_multi_bb_lsm_multidim.stan")
i.lsm.mbb.fit <- int_multi_bb_lsm$sample(
  data = dat.list.int, seed = 12345, iter_warmup = 750,
  iter_sampling = 750, chains = 3, parallel_chains = 3,
  init = function () list(
    lambda_p = rep(.5, dat.list.int$Ni), p_lgt = rep(0, dat.list.int$Ni),
    r_lgt = rep(0, dat.list.int$Ni), sigma_eta = .5))
i.lsm.mbb.fit$cmdstan_diagnose()
i.lsm.mbb.fit <- read_stan_csv(i.lsm.mbb.fit$output_files())
print(i.lsm.mbb.fit, c("Eta_pr"), digits_summary = 3, include = FALSE)
# saveRDS(i.lsm.mbb.fit, "sim_dat_res/i.lsm.mbb.fit.md.rds")
# i.lsm.mbb.fit <- readRDS("sim_dat_res/i.lsm.mbb.fit.md.rds")
rowSums(get_elapsed_time(i.lsm.mbb.fit))
# chain:1 chain:2 chain:3 
# 869.462 842.473 862.813 

print(mbb.lsm.fit, c("R[1,2]", "R[1,3]", "R[2,3]", "sigma_eta"), digits_summary = 3,
      probs = c(.025, .05, .25, .5, .75, .95, .975))
#             mean se_mean    sd   2.5%     5%    25%    50%    75%    95%  97.5% n_eff  Rhat
# R[1,2]     0.521   0.003 0.045  0.429  0.445  0.493  0.522  0.553  0.591  0.603   171 1.011
# R[1,3]    -0.291   0.006 0.069 -0.428 -0.406 -0.335 -0.291 -0.244 -0.177 -0.155   132 1.015
# R[2,3]    -0.244   0.005 0.068 -0.373 -0.352 -0.291 -0.248 -0.198 -0.128 -0.105   177 1.027
# sigma_eta  1.062   0.003 0.066  0.938  0.952  1.018  1.062  1.106  1.172  1.191   649 1.007

print(i.lsm.mbb.fit, c("R[1,2]", "R[1,3]", "R[2,3]", "sigma_eta"), digits_summary = 3,
      probs = c(.025, .05, .25, .5, .75, .95, .975))
#             mean se_mean    sd   2.5%     5%    25%    50%    75%    95%  97.5% n_eff  Rhat
# R[1,2]     0.523   0.002 0.047  0.429  0.445  0.492  0.524  0.556  0.597  0.613   427 1.003
# R[1,3]    -0.298   0.004 0.085 -0.477 -0.446 -0.352 -0.297 -0.240 -0.164 -0.138   428 1.010
# R[2,3]    -0.218   0.004 0.086 -0.387 -0.360 -0.274 -0.216 -0.161 -0.084 -0.058   369 1.006
# sigma_eta  0.984   0.005 0.093  0.811  0.838  0.922  0.982  1.044  1.138  1.167   372 1.011

params$means
print(i.lsm.mbb.fit, c("p"), digits_summary = 3)

params$iccs
print(i.lsm.mbb.fit, c("rho"), digits_summary = 3)

params$lambda
print(i.lsm.mbb.fit, L.idxs, digits_summary = 3)

lsm.fit.df <- as.data.frame(summary(
  mbb.lsm.fit, c("R[1,2]", "R[1,3]", "R[2,3]", "sigma_eta", L.idxs,
                 "p", "rho"), probs = c(.05, .5, .95))$summary)
lsm.fit.df$params <- rownames(lsm.fit.df)
i.lsm.fit.df <- as.data.frame(summary(
  i.lsm.mbb.fit, c("R[1,2]", "R[1,3]", "R[2,3]", "sigma_eta", L.idxs,
                   "p", "rho"), probs = c(.05, .5, .95))$summary)
i.lsm.fit.df$params <- rownames(i.lsm.fit.df)

comb.df <- rbindlist(list(lsm.fit.df, i.lsm.fit.df), idcol = "fit")
comb.df
comb.df[, class := gsub("\\[\\d(,\\d)?\\]", "", params)]
comb.df[, class := ifelse(class == "Lambda_mat", "loading", class)]
comb.df[, class := ifelse(class == "p", "mean", class)]
comb.df[, class := ifelse(class == "rho", "ICC", class)]
comb.df[, class := ifelse(class == "sigma_eta", "load_icc", class)]
comb.df[, item := as.integer(gsub("[a-z]+|[A-Z]+|\\[|\\]|\\_|,\\d", "", params))]
setnames(comb.df, c("50%", "5%", "95%"), c("median", "ll", "ul"))
comb.df

comb.df[, pos := ifelse(
  class == "loading", as.integer(item), ifelse(
    class == "mean", as.integer(item) + 6, ifelse(
      class == "ICC", as.integer(item) + 12, 13)))]

params$item <- params$time

comb.df <- merge(comb.df, params, all.x = TRUE)
comb.df

comb.df[, pop_param := ifelse(
  class == "mean", means, ifelse(
    class == "loading", abs(lambdas), ifelse(
      class == "ICC", iccs, ifelse(
        class == "load_ICC", 1, NA))))]

comb.df[, fit.t := factor(
  fit, labels = c("Location-ICC", "Location-ICC (binned data)"))]
comb.df

comb.df[pop_param == .9 & class == "mean", median := 1 - median]
comb.df[pop_param == .9 & class == "mean", ll := 1 - ll]
comb.df[pop_param == .9 & class == "mean", ul := 1 - ul]
comb.df[pop_param == .9 & class == "mean", pop_param := 1 - pop_param]
comb.df

library(ggplot2)
library(ggforce)
library(scales)

# Fits: 1 = Location-ICC, 2 = Location-ICC (binned data)

# More extreme loadings from location-only fit
# Smaller ICCs (in high ICC models) from location-only fit
plt.1 <- ggplot(comb.df[class %in% c("loading", "mean", "ICC")],
       aes(item, abs(median), fill = fit.t)) +
  geom_segment(aes(x = item - .25, xend = item + .25, y = pop_param, yend = pop_param)) +
  geom_linerange(aes(ymin = abs(ll), ymax = abs(ul), col = fit.t),
                 position = position_dodge(.5)) +
  facet_wrap(~ reorder(class, pos), scales = "free") +
  theme_bw() +
  scale_x_continuous(labels = 1:7, breaks = 1:7) +
  scale_color_manual(values = cbPalette) +
  theme(legend.position = "top", strip.background = element_blank(),
        panel.border = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), panel.grid.major.x = element_blank()) +
  labs(x = "Item ID", y = "90% quantile intervals", fill = "", col = "")

cbrdf <- comb.df[class == "R"]
cbrdf$class <- "|R|"
cbrdf$pop_param <- c(rep(c(.5, .35), 2), rep(.25, 2))
cbrdf$item <- c(rep(c(1, 2), 2), rep(3, 2))
cbrdf$item.f <- c(rep(c("F1-F2", "F1-F3"), 2), rep("F2-F3", 2))

plt.2 <- ggplot(cbrdf, aes(item, abs(median), fill = fit.t)) +
  geom_linerange(aes(ymin = abs(ll), ymax = abs(ul), col = fit.t),
                 position = position_dodge(.5)) +
  geom_segment(aes(x = item - .25, xend = item + .25, y = pop_param, yend = pop_param)) +
  facet_wrap(~ reorder(class, pos), scales = "free") +
  scale_x_continuous(labels = c("F1-F2", "F1-F3", "F2-F3"), breaks = 1:3) +
  scale_y_continuous(breaks = c(.1, .25, .35, .5, .6)) +
  theme_bw() + guides(col = FALSE) +
  scale_color_manual(values = cbPalette) +
  theme(legend.position = "top", strip.background = element_blank(),
        panel.border = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), panel.grid.major.x = element_blank(),
        axis.title.y = element_blank()) +
  labs(x = "Correlation", y = "90% quantile intervals", fill = "", col = "")

library(patchwork)

plt.1 + plt.2 + plot_layout(widths = c(3, 1))
ggsave("plots/sim_2_res.pdf", height = 3.5, width = 6.5)
ggsave(paste0(plot_folder, "sim_2_res.pdf"), height = 3.5, width = 6.5)

