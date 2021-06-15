# SIMULATION DATA ANALYSIS # 1

plot_folder <- "~/Dropbox/dissertation/paper_3/plots/"

library(MASS)
library(rstan)
library(cmdstanr)
library(data.table)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

set_cmdstan_path("~/cmdstan/")

beta_bin <- cmdstan_model("stan_scripts/multi_beta_bin.stan")
beta_bin_lsm <- cmdstan_model("stan_scripts/multi_beta_bin_lsm.stan")
beta_bin_interval <- cmdstan_model("stan_scripts/interval_multi_bb.stan")
beta_bin_lsm_interval <- cmdstan_model("stan_scripts/interval_multi_bb_lsm.stan")

(params <- expand.grid(
  lambdas = c(.35, 1), means = c(.05, .2, .9), iccs = c(.8, .5, .2)))
params$precs <- exp(-qlogis(params$iccs))
params$lambdas[params$means == .9] <- params$lambdas[params$means == .9] * -1
params

n.i <- nrow(params)
r.traits <- -.35

set.seed(12345)
cov(Eta <- mvrnorm(2e3, c(0, 0), matrix(c(1, -r.traits, -r.traits, 1), 2), empirical = TRUE))
plot(Eta)
cov(M1 <- t(matrix(params$lambda) %*% Eta[, 1] + qlogis(params$means)))
plogis(colMeans(M1))
cov(M2 <- t(matrix(rep(1, n.i)) %*% Eta[, 2] + log(params$precs)))
exp(colMeans(M2))

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

# ICC/prec not count for these plots
{
  png(paste0(plot_folder, "sim_prob_dist.png"), width = 6.5, height = 3.5, units = "in", res = 360)
  par(mfrow = c(2, 3))
  for (j in c(1, 3, 5, 2, 4, 6)) hist(
    plogis(M1[, j]), main = paste0("m=",params$means[j], ",l=", params$lambdas[j]),
    xlab = "probability")
  rm(j)
  par(mfrow = c(1, 1))
  dev.off()
}

# All parameters count for these plots
par(mfrow = c(2, 4))
for (j in 1:ncol(N)) {
  if (params$iccs[j] == .5 | params$means[j] == .9) next
  barplot(table(N[, j]), main = paste0("i=",params$iccs[j], ",m=",params$means[j], ",l=", params$lambdas[j]),
          xlab = "count", ylim = c(0, 1700))
}
rm(j)
par(mfrow = c(1, 1))

# Not so important plots
par(mfrow = c(3, 3))
for (j in 1:ncol(M2)) hist(exp(M2[, j]), main = params$precs[j])
rm(j)
par(mfrow = c(1, 1))

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
N.df.dt

library(scales)

ggplot(N.df.dt[iccs != .5], aes(n, N)) +
  geom_bar(aes(fill = factor(abs(lambdas))), col = 1, stat = "identity", position = position_dodge(.5)) +
  facet_grid(iccs.t ~ means.t) +
  theme_bw() + theme(legend.position = "top") +
  labs(x = "Counts", y = "Frequency", fill = "Absolute value of loading") +
  # scale_fill_manual(values = cbPalette) +
  scale_fill_grey(start = .8, end = 0)
ggsave(paste0(plot_folder, "sim_item_dist.pdf"), width = 6.5, height = 4.5)

dat.list <- list(
  Np = nrow(N), Ni = ncol(N), NpNi = nrow(N.df), N = 7,
  n_1d = N.df$n, p_ids = N.df$id, i_ids = N.df$time, shape_r = 2,
  lambda_median = log(.5), lambda_scale = log(4 / .5) / qnorm(.99),
  n_neg = sum(params$lambda < 0), neg_l = which(params$lambda < 0))
# saveRDS(dat.list, "sim_dat_res/sim_dat_list.rds")
# dat.list <- readRDS("sim_dat_res/sim_dat_list.rds")

multi_beta_bin <- cmdstan_model("stan_scripts/multi_beta_bin.stan")
mbb.fit <- multi_beta_bin$sample(
  data = dat.list, seed = 12345, iter_warmup = 5e2,
  iter_sampling = 5e2, chains = 3, parallel_chains = 3)
mbb.fit$cmdstan_diagnose()
mbb.fit <- read_stan_csv(mbb.fit$output_files())
print(mbb.fit, c("eta"), digits_summary = 3, include = FALSE)
# saveRDS(mbb.fit, "sim_dat_res/mbb.fit.rds")
# mbb.fit <- readRDS("sim_dat_res/mbb.fit.rds")
rowSums(get_elapsed_time(mbb.fit))
# chain:1 chain:2 chain:3
# 354.373 271.996 356.379

params$means
print(mbb.fit, c("p"), digits_summary = 3)

params$iccs
print(mbb.fit, c("rho"), digits_summary = 3)

params$lambda
print(mbb.fit, c("lambda"), digits_summary = 3)

multi_beta_bin_lsm <- cmdstan_model("stan_scripts/multi_beta_bin_lsm.stan")
mbb.fit.lsm <- multi_beta_bin_lsm$sample(
  data = dat.list, seed = 123456, iter_warmup = 5e2,
  iter_sampling = 5e2, chains = 3, parallel_chains = 3)
mbb.fit.lsm$cmdstan_diagnose()
mbb.fit.lsm <- read_stan_csv(mbb.fit.lsm$output_files())
print(mbb.fit.lsm, c("Eta"), digits_summary = 3, include = FALSE)
# saveRDS(mbb.fit.lsm, "sim_dat_res/mbb.fit.lsm.rds")
# mbb.fit.lsm <- readRDS("sim_dat_res/mbb.fit.lsm.rds")
rowSums(get_elapsed_time(mbb.fit.lsm))
# chain:1 chain:2 chain:3
# 588.810 589.594 612.272

print(mbb.fit.lsm, c("r", "sigma_eta"), digits_summary = 3,
      probs = c(.025, .05, .25, .5, .75, .95, .975))
#             mean se_mean    sd   2.5%     5%    25%    50%    75%    95%  97.5% n_eff  Rhat
# r         -0.373   0.002 0.035 -0.440 -0.431 -0.398 -0.374 -0.350 -0.314 -0.303   441 1.004
# sigma_eta  1.038   0.001 0.033  0.978  0.984  1.016  1.039  1.061  1.091  1.103   591 1.007

params$means
print(mbb.fit.lsm, c("p"), digits_summary = 3)

params$iccs
print(mbb.fit.lsm, c("rho"), digits_summary = 3)

params$lambda
print(mbb.fit.lsm, c("lambda"), digits_summary = 3)

table(N.df$q <- cut(dat.list$n_1d, c(-1, 0, 2, 4, 7)))
table(N.df$q.int <- as.integer(N.df$q))

N.df.dt <- as.data.table(N.df)[, .N, list(time, q.int)]
N.df.dt

dat.list.mini <- list(
  N = 7, n_1d = N.df.dt$q.int, i_ids = N.df.dt$time, counts = N.df.dt$N)
dat.list.mini$Ni <- max(dat.list.mini$i_ids)
dat.list.mini$Nc <- length(dat.list.mini$n_1d)
dat.list.mini$N_int <- 4
# first cut system
(dat.list.mini$cuts <- matrix(c(0, 0, 1, 2, 3, 4, 5, 7), byrow = TRUE, ncol = 2))

int_bb <- cmdstan_model("stan_scripts/interval_bb_item.stan")
i.bb.fit <- int_bb$optimize(data = dat.list.mini)
(ll.bb <- sum(i.bb.fit$summary()[regexpr("log_lik", i.bb.fit$summary()$variable) > 0, 2]))

# second cut system
(dat.list.mini$cuts <- matrix(c(0, 1, 1, 2, 3, 4, 5, 7), byrow = TRUE, ncol = 2))
i.bb.fit.2 <- int_bb$optimize(data = dat.list.mini)
(ll.bb.2 <- sum(i.bb.fit.2$summary()[regexpr("log_lik", i.bb.fit.2$summary()$variable) > 0, 2]))

# ordinal system
int_ord <- cmdstan_model("stan_scripts/ordinal_item.stan")
i.ord.fit <- int_ord$optimize(data = dat.list.mini)
(ll.ord <- sum(i.ord.fit$summary()[regexpr("log_lik", i.ord.fit$summary()$variable) > 0, 2]))

round(c(bb = 2 * (dat.list.mini$Ni * 2) - 2 * ll.bb,
        bb.2 = 2 * (dat.list.mini$Ni * 2) - 2 * ll.bb.2,
        ord = 2 * (dat.list.mini$Ni * 3) - 2 * ll.ord))

dat.list.int <- dat.list
dat.list.int$n_1d <- N.df$q.int
dat.list.int$N_int <- 4
# dat.list.int$cuts <- c(0, 0, 2, 4, 7)
(dat.list.int$cuts <- matrix(c(0, 0, 1, 2, 3, 4, 5, 7), byrow = TRUE, ncol = 2))
dat.list.int$count_1 <- sum(dat.list.int$n_1d == 1)
dat.list.int$count_max <- sum(dat.list.int$n_1d == 4)
dat.list.int$pos_1 <- which(dat.list.int$n_1d == 1)
dat.list.int$pos_max <- which(dat.list.int$n_1d == 4)

dat.list.int$cuts
table(dat.list$n_1d, dat.list.int$n_1d)

int_multi_bb <- cmdstan_model("stan_scripts/interval_multi_bb.stan")
i.mbb.fit <- int_multi_bb$sample(
  data = dat.list.int, seed = 123, iter_warmup = 5e2,
  iter_sampling = 5e2, chains = 3, parallel_chains = 3,
  init = function () list(
    lambda_p = rep(.5, dat.list.int$Ni), p_lgt = rep(0, dat.list.int$Ni),
    r_lgt = rep(0, dat.list.int$Ni), sigma_eta = .5))
i.mbb.fit$cmdstan_diagnose()
i.mbb.fit <- read_stan_csv(i.mbb.fit$output_files())
print(i.mbb.fit, c("eta"), digits_summary = 3, include = FALSE)
# saveRDS(i.mbb.fit, "sim_dat_res/i.mbb.fit.rds")
# i.mbb.fit <- readRDS("sim_dat_res/i.mbb.fit.rds")
rowSums(get_elapsed_time(i.mbb.fit))
#  chain:1  chain:2  chain:3 
# 1113.872 1082.994 1121.843 

params$means
print(i.mbb.fit, c("p"), digits_summary = 3)

params$iccs
print(i.mbb.fit, c("rho"), digits_summary = 3)

params$lambda
print(i.mbb.fit, c("lambda"), digits_summary = 3)

int_multi_bb_lsm <- cmdstan_model("stan_scripts/interval_multi_bb_lsm.stan")
i.lsm.mbb.fit <- int_multi_bb_lsm$sample(
  data = dat.list.int, seed = 1, iter_warmup = 5e2,
  iter_sampling = 5e2, chains = 3, parallel_chains = 3,
  init = function () list(
    lambda_p = rep(.5, dat.list.int$Ni), p_lgt = rep(0, dat.list.int$Ni),
    r_lgt = rep(0, dat.list.int$Ni), sigma_eta = .5, r = 0))
i.lsm.mbb.fit$cmdstan_diagnose()
i.lsm.mbb.fit <- read_stan_csv(i.lsm.mbb.fit$output_files())
print(i.lsm.mbb.fit, c("Eta"), digits_summary = 3, include = FALSE)
# saveRDS(i.lsm.mbb.fit, "sim_dat_res/i.lsm.mbb.fit.rds")
# i.lsm.mbb.fit <- readRDS("sim_dat_res/i.lsm.mbb.fit.rds")
rowSums(get_elapsed_time(i.lsm.mbb.fit))
#  chain:1  chain:2  chain:3 
# 1432.267 1477.254 1452.091 

print(mbb.fit.lsm, c("r", "sigma_eta"), digits_summary = 3,
      probs = c(.025, .05, .25, .5, .75, .95, .975))
#             mean se_mean    sd   2.5%     5%    25%    50%    75%    95%  97.5% n_eff  Rhat
# r         -0.373   0.002 0.035 -0.440 -0.431 -0.398 -0.374 -0.350 -0.314 -0.303   441 1.004
# sigma_eta  1.038   0.001 0.033  0.978  0.984  1.016  1.039  1.061  1.091  1.103   591 1.007

print(i.lsm.mbb.fit, c("r", "sigma_eta"), digits_summary = 3,
      probs = c(.025, .05, .25, .5, .75, .95, .975))
#             mean se_mean    sd   2.5%     5%    25%    50%    75%    95%  97.5% n_eff  Rhat
# r         -0.383   0.002 0.041 -0.462 -0.447 -0.412 -0.384 -0.355 -0.315 -0.301   491 1.007
# sigma_eta  1.060   0.002 0.039  0.981  0.995  1.034  1.062  1.086  1.122  1.133   639 1.001

params$means
print(i.lsm.mbb.fit, c("p"), digits_summary = 3)

params$iccs
print(i.lsm.mbb.fit, c("rho"), digits_summary = 3)

params$lambda
print(i.lsm.mbb.fit, c("lambda"), digits_summary = 3)

base.fit.df <- as.data.frame(summary(
  mbb.fit, c("lambda", "p", "rho"), probs = c(.05, .5, .95))$summary)
base.fit.df$params <- rownames(base.fit.df)
lsm.fit.df <- as.data.frame(summary(
  mbb.fit.lsm, c("r", "sigma_eta", "lambda", "p", "rho"), probs = c(.05, .5, .95))$summary)
lsm.fit.df$params <- rownames(lsm.fit.df)
i.fit.df <- as.data.frame(summary(
  i.mbb.fit, c("lambda", "p", "rho"), probs = c(.05, .5, .95))$summary)
i.fit.df$params <- rownames(i.fit.df)
i.lsm.fit.df <- as.data.frame(summary(
  i.lsm.mbb.fit, c("r", "sigma_eta", "lambda", "p", "rho"), probs = c(.05, .5, .95))$summary)
i.lsm.fit.df$params <- rownames(i.lsm.fit.df)

comb.df <- rbindlist(list(base.fit.df, lsm.fit.df, i.fit.df, i.lsm.fit.df), idcol = "fit")
comb.df
comb.df[, class := gsub("\\[\\d+\\]", "", params)]
comb.df[, class := ifelse(class == "lambda", "loading", class)]
comb.df[, class := ifelse(class == "p", "mean", class)]
comb.df[, class := ifelse(class == "rho", "ICC", class)]
comb.df[, class := ifelse(class == "sigma_eta", "load_icc", class)]
comb.df[, item := as.integer(gsub("[a-z]+|\\[|\\]|\\_", "", params))]
setnames(comb.df, c("50%", "5%", "95%"), c("median", "ll", "ul"))
comb.df

comb.df[, params := ifelse(is.na(item), class, paste0(class, "[", item, "]"))]
comb.df
comb.df[, pos := ifelse(
  class == "loading", as.integer(item), ifelse(
    class == "mean", as.integer(item) + 18, ifelse(
      class == "ICC", as.integer(item) + 36, 55)))]

comb.df[, params := ifelse(item %in% dat.list$neg_l & class == "loading", paste0("-", params), params)]

params$item <- params$time

comb.df <- merge(comb.df, params, all.x = TRUE)
comb.df

comb.df[, pop_param := ifelse(
  class == "mean", means, ifelse(
    class == "loading", abs(lambdas), ifelse(
      class == "ICC", iccs, ifelse(
        class == "load_ICC", 1, .35))))]

comb.df[, grouper := paste0(class, " = ", pop_param)]
comb.df

comb.df[, fit.t := factor(
  fit, labels = c("Location", "Location-ICC", "Location (binned data)", "Location-ICC (binned data)"))]
comb.df

library(ggplot2)
library(ggforce)
library(scales)

# Fits: 1 = Location, 2 = Location-ICC, 3 = Location (binned data), 4 = Location-ICC (binned data)

# More extreme loadings from location-only fit
# Smaller ICCs (in high ICC models) from location-only fit
ggplot(comb.df[!is.na(item)],
       aes(reorder(item, pos), abs(median), fill = fit.t)) +
  geom_hline(aes(yintercept = pop_param), size = .5, linetype = 1, alpha = .5) +
  geom_linerange(aes(ymin = abs(ll), ymax = abs(ul), col = fit.t),
                 position = position_dodge(.5)) +
  facet_wrap(~ reorder(grouper, pos), scales = "free", strip.position = "top", ncol = 2) +
  theme_bw() +
  scale_color_manual(values = cbPalette) +
  # scale_color_grey(start = .8, end = 0) +
  scale_shape_manual(values = c(1, 5, 4, 8)) +
  theme(legend.position = "top", strip.background = element_blank(),
        panel.border = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), panel.grid.major.x = element_blank()) +
  labs(x = "Item ID", y = "90% quantile intervals", fill = "", col = "",
       caption = "Values for loadings 3, 6 & 9 are sign reversed")
ggsave("plots/sim_1_res.pdf", height = 6, width = 6.5)
ggsave(paste0(plot_folder, "sim_1_res.pdf"), height = 6, width = 6.5)

