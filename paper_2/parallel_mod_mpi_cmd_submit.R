data.location <- "summary_text/mpi_dat/"  # location of csv dataset

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Load packages (Check for packages that should be installed) ----
source("parallel_mod_mpi_cmd_load_packages.R")

# Pre-process data (Worth reviewing source file) ----
source("parallel_mod_mpi_cmd_prep_data.R")

stan.scripts <- "summary_text/stan_scripts/"  # location of stan scripts
print.images <- "summary_text/text_clean/"  # location to print images

# FIGURE 1 ----
ggplot(X, aes(response.l)) + geom_histogram(col = 1) + facet_wrap(~ item_n) +
  scale_x_continuous(name = "Probit-transformed prevalence of deprivation") +
  theme(strip.background = element_blank(), panel.border = element_blank(),
        panel.spacing.x = unit(.5, "cm"), panel.grid.minor.y = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = margin(0, 0, 0, 0)) +
  labs(y = "Number of countries")
ggsave(paste0(print.images, "01_deprivation.pdf"), width = 6, height = 2.75)

# Table 1 ----
describe(dat[, 15:20])
#    vars   n  mean   sd median trimmed  mad   min  max range  skew kurtosis   se
# x1    1 100 -0.06 1.45  -0.17   -0.11 1.64 -2.34 2.58  4.91  0.17    -1.01 0.14
# x2    2 101  0.31 1.03   0.21    0.30 1.20 -1.95 2.24  4.19  0.10    -1.01 0.10
# x3    3 101  0.80 0.77   0.67    0.77 0.78 -0.53 2.72  3.26  0.36    -0.83 0.08
# x4    4 100  0.85 1.29   0.86    0.90 1.86 -1.66 2.58  4.24 -0.16    -1.37 0.13
# x5    5 100  0.38 1.04   0.33    0.40 1.21 -2.10 2.40  4.50 -0.19    -0.65 0.10
# x6    6 101  1.31 0.96   1.13    1.29 1.09 -0.53 3.50  4.02  0.17    -0.89 0.10

round(cor(dat[, 15:20], use = "p"), 2)

# MODELLING ----

# CONGENERIC CFA ----
# Congeneric CFA in lavaan ----
# Creating error SD as square root of error variance
X.tmp <- dat[, 15:20]
colnames(X.tmp) <- c("cf", "sanit", "dw", "elect", "hous", "ast")
writeLines(lav.form <- paste(
  "F =~ cf + sanit + dw + elect + hous + ast
  cf ~~ a * cf
  sanit ~~ b * sanit
  dw ~~ c * dw
  elect ~~ d * elect
  hous ~~ e * hous
  ast ~~ f * ast",
  "a1 := sqrt(a)\nb1 := sqrt(b)\nc1 := sqrt(c)\nd1 := sqrt(d)\ne1 := sqrt(e)\nf1 := sqrt(f)",
  sep = "\n"))
summary(lav.fit <- cfa(
  lav.form, X.tmp, std.lv = TRUE, meanstructure = TRUE, missing = "ML"),
  fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
rm(X.tmp)

# FIGURE 2 ----
# plot fit
png(paste0(print.images, "02_lavfit.png"), width = 5, height = 2, units = "in", res = 150)
semPaths(lav.fit, "est", fade = FALSE, style = "ram", edge.label.cex = 1.5,
         mar = c(2, 2, 2, 2), edge.color = cbPalette[1])
dev.off()

parameterEstimates(lav.fit, standardized = TRUE, rsquare = TRUE)
#      lhs op     rhs label    est    se       z pvalue ci.lower ci.upper std.lv std.all std.nox
# 1      F =~      cf        1.280 0.113  11.361  0.000    1.059    1.501  1.280   0.888   0.888
# 2      F =~   sanit        0.926 0.078  11.791  0.000    0.772    1.079  0.926   0.906   0.906
# 3      F =~      dw        0.623 0.063   9.897  0.000    0.499    0.746  0.623   0.812   0.812
# 4      F =~   elect        1.248 0.093  13.453  0.000    1.066    1.429  1.248   0.974   0.974
# 5      F =~    hous        0.883 0.083  10.684  0.000    0.721    1.044  0.883   0.855   0.855
# 6      F =~     ast        0.882 0.073  12.142  0.000    0.740    1.025  0.882   0.921   0.921
# 14    cf ~1               -0.043 0.144 -0.296  0.767   -0.324    0.239 -0.043  -0.029  -0.029
# 15 sanit ~1                0.313 0.102  3.078  0.002    0.114    0.512  0.313   0.306   0.306
# 16    dw ~1                0.800 0.076 10.490  0.000    0.651    0.949  0.800   1.044   1.044
# 17 elect ~1                0.849 0.127  6.661  0.000    0.599    1.099  0.849   0.663   0.663
# 18  hous ~1                0.388 0.103  3.767  0.000    0.186    0.589  0.388   0.375   0.375
# 19   ast ~1                1.307 0.095 13.710  0.000    1.120    1.493  1.307   1.364   1.364
# 21    a1 := sqrt(a)    a1  0.664 0.052  12.700  0.000    0.562    0.766  0.664   0.460   0.460
# 22    b1 := sqrt(b)    b1  0.433 0.035  12.314  0.000    0.364    0.502  0.433   0.424   0.424
# 23    c1 := sqrt(c)    c1  0.447 0.033  13.428  0.000    0.382    0.512  0.447   0.583   0.583
# 24    d1 := sqrt(d)    d1  0.288 0.043   6.753  0.000    0.204    0.372  0.288   0.225   0.225
# 25    e1 := sqrt(e)    e1  0.536 0.041  13.048  0.000    0.456    0.617  0.536   0.519   0.519
# 26    f1 := sqrt(f)    f1  0.372 0.031  11.851  0.000    0.311    0.434  0.372   0.389   0.389
# 27    cf r2      cf        0.788    NA      NA     NA       NA       NA     NA      NA      NA
# 28 sanit r2   sanit        0.821    NA      NA     NA       NA       NA     NA      NA      NA
# 29    dw r2      dw        0.660    NA      NA     NA       NA       NA     NA      NA      NA
# 30 elect r2   elect        0.949    NA      NA     NA       NA       NA     NA      NA      NA
# 31  hous r2    hous        0.730    NA      NA     NA       NA       NA     NA      NA      NA
# 32   ast r2     ast        0.849    NA      NA     NA       NA       NA     NA      NA      NA

# FIGURE 3 ----
# plot log-normal
(loc.ln <- log(sqrt(.3)))
(scale.ln <- log(sqrt(3) / sqrt(.3)) / qnorm(.95))
(ggplot(data.frame(x = seq(0, 2, 1e-4), density = dlnorm(seq(0, 2, 1e-4), loc.ln, scale.ln)),
        aes(x, density)) + geom_line() + theme_classic() +
   scale_x_continuous(breaks = c(0, exp(loc.ln), sqrt(3)), labels = number_format()) +
   geom_vline(xintercept = c(exp(loc.ln), sqrt(3)), linetype = 2, alpha = .5) +
   theme(axis.title = element_blank()) +
   labs(tag = "A", subtitle = "Density")) +
  (ggplot(data.frame(x = seq(0, 2, 1e-4), cum_dens = plnorm(seq(0, 2, 1e-4), loc.ln, scale.ln)),
          aes(x, cum_dens)) + geom_line() + theme_classic() +
     scale_x_continuous(breaks = c(0, exp(loc.ln), sqrt(3)), labels = number_format()) +
     scale_y_continuous(labels = percent_format(), breaks = c(.0, .25, .5, .75, .95)) +
     theme(axis.title.x = element_text(hjust = 1.5), axis.title.y = element_blank()) +
     geom_hline(yintercept = c(.5, .95), linetype = 2, alpha = .5) +
     geom_vline(xintercept = c(exp(loc.ln), sqrt(3)), linetype = 2, alpha = .5) +
     labs(x = "A-priori expectations for lambda_i parameter",
          subtitle = "Cumulative distribution", tag = "B"))
ggsave(paste0(print.images, "03_log_normal_lambda.pdf"), width = 5, height = 2.5)

# Congeneric regression in Stan ----
# Data list with references to equation 11
data <- list(
  alpha_scale = 3,  # scale of a parameter
  beta_scale = 2,  # scale of sigma_t
  sigma_scale = sqrt(3),  # scale of sigma_i
  lambda_median = loc.ln,  # log-median of lambda_i
  lambda_scale = scale.ln,  # scale of lambda_i
  # item id number long-form, ids be contiguous from 1 to number of items
  item_id = X$item,
  Ni = max(X$item),  # number of items
  # respondent id number long-form, ids be contiguous from 1 to number of respondents
  resp_id = X$id,
  Np = max(X$id),  # number of respondents
  y = X$response.l,  # logit-response
  N = length(X$response),  # rows in data (long form)
  ret_yhat = 1,  # return predicted outcomes (1: yes, 0: no)
  ret_ll = 0  # return log-likelihood (1: yes, 0 = no),
)

# Compile Stan model (takes time)
mod <- cmdstan_model(file.path(paste0(stan.scripts, "cong.stan")))
fit.0 <- mod$sample(
  data = data,  # data list passed to Stan
  seed = 12345,  # Seed for reproducibility
  iter_warmup = 1e3,  # 1000 warmup samples per chain
  iter_sampling = 3e3,  # 3000 samples for inference per chain
  chains = 8,  # 8 chains
  parallel_chains = 8  # 8 cores on multicore systems,
)
fit.0$cmdstan_diagnose()  # Run quick diagnostic checks
fit.0 <- read_stan_csv(fit.0$output_files())  # convert to Rstan

# Print major parameters
print(fit.0, c("theta_p", "yhat", "lp__"), include = FALSE)
#              mean se_mean   sd  2.5%   25%   50%   75% 97.5% n_eff Rhat
# alpha        0.55    0.02 0.79 -1.03  0.01  0.56  1.09  2.09  2647 1.00
# beta_dev[1] -0.59    0.02 0.79 -2.13 -1.12 -0.59 -0.06  0.97  2669 1.00
# beta_dev[2] -0.23    0.02 0.79 -1.77 -0.77 -0.24  0.30  1.33  2674 1.00
# beta_dev[3]  0.25    0.02 0.79 -1.27 -0.28  0.25  0.78  1.81  2680 1.00
# beta_dev[4]  0.31    0.02 0.79 -1.23 -0.23  0.30  0.83  1.86  2657 1.00
# beta_dev[5] -0.16    0.02 0.79 -1.69 -0.70 -0.17  0.37  1.41  2677 1.00
# beta_dev[6]  0.76    0.02 0.79 -0.77  0.22  0.75  1.29  2.32  2673 1.00
# lambda[1]    1.24    0.00 0.11  1.04  1.17  1.24  1.32  1.47  2157 1.00
# lambda[2]    0.90    0.00 0.08  0.77  0.85  0.90  0.95  1.06  2065 1.00
# lambda[3]    0.61    0.00 0.06  0.49  0.56  0.60  0.65  0.73  2729 1.00
# lambda[4]    1.22    0.00 0.09  1.06  1.16  1.22  1.28  1.40  1517 1.00
# lambda[5]    0.86    0.00 0.08  0.71  0.80  0.86  0.91  1.03  2428 1.00
# lambda[6]    0.86    0.00 0.07  0.73  0.81  0.86  0.91  1.01  1891 1.00
# sigma[1]     0.68    0.00 0.05  0.58  0.64  0.68  0.72  0.79 16998 1.00
# sigma[2]     0.44    0.00 0.04  0.38  0.42  0.44  0.47  0.52 16476 1.00
# sigma[3]     0.46    0.00 0.04  0.39  0.43  0.46  0.48  0.53 21096 1.00
# sigma[4]     0.29    0.00 0.05  0.20  0.26  0.29  0.32  0.38  3977 1.00
# sigma[5]     0.55    0.00 0.04  0.47  0.52  0.55  0.58  0.64 18874 1.00
# sigma[6]     0.38    0.00 0.03  0.32  0.36  0.38  0.40  0.45 14102 1.00
# beta[1]     -0.04    0.00 0.14 -0.31 -0.13 -0.04  0.06  0.24  1262 1.00
# beta[2]      0.32    0.00 0.10  0.12  0.25  0.32  0.38  0.51  1218 1.01
# beta[3]      0.80    0.00 0.08  0.66  0.75  0.80  0.85  0.95  1520 1.00
# beta[4]      0.86    0.00 0.12  0.61  0.77  0.85  0.94  1.10  1021 1.01
# beta[5]      0.39    0.00 0.10  0.19  0.32  0.39  0.46  0.59  1350 1.01
# beta[6]      1.31    0.00 0.09  1.13  1.25  1.31  1.37  1.50  1164 1.01

# saveRDS(fit.0, "res/mod_1.rds")
fit.0 <- readRDS("res/mod_1.rds")
check_hmc_diagnostics(fit.0)
# 2 of 24000 iterations ended with a divergence (0.00833333333333333%).

# coefficient omega:
sum(apply(as.data.frame(fit.0, "lambda"), 2, median)) ^ 2 / (
  sum(apply(as.data.frame(fit.0, "lambda"), 2, median)) ^ 2 +
    sum(apply(as.data.frame(fit.0, "sigma"), 2, median) ^ 2)
)

# Traceplots:
mcmc_trace(fit.0, pars = c(
  paste0("beta[", 1:6, "]"), paste0("lambda[", 1:6, "]"),
  paste0("sigma[", 1:6, "]")))
# Rank plots:
mcmc_rank_overlay(
  fit.0, pars = c(
    paste0("beta[", 1:6, "]"), paste0("lambda[", 1:6, "]"),
    paste0("sigma[", 1:6, "]")), ref_line = TRUE) +
  scale_x_continuous(labels = number_format(scale = 1e-3, 1), name = "Rank (1000s)") +
  theme(legend.position = "top", strip.background = element_blank())
ggsave(paste0(print.images, "app_01_rank_mod_1.pdf"), width = 6.5, height = 5)

pars.0.df <- as.data.frame(
  summary(fit.0, c("beta", "lambda", "sigma"))$summary)
# beta[1-6] are means
# lambda[1-6] are loadings
# sigma[1-6] are residual standard deviations
pars.0.df$variable <- rownames(pars.0.df)
# Join lavaan and Bayes results
pars.0.df
pars.0.df.joined <- as.data.table(cbind(
  pars.0.df[, c(11, 6, 4, 8)],
  parameterEstimates(lav.fit)[c(14:19, 1:6, 21:26), c(5, 9:10)]))
pars.0.df.joined
# Rename columns
setnames(pars.0.df.joined, 2:7, c("m.1", "ll.1", "ul.1", "m.2", "ll.2", "ul.2"))
pars.0.df.joined <- reshape(
  pars.0.df.joined, direction = "long", varying = list(c(2, 5), c(3, 6), c(4, 7)))
pars.0.df.joined[, model := c("Bayesian", "ML")[time]]
pars.0.df.joined[, type := factor(
  gsub("\\[\\d\\]", "", variable), levels = c("lambda", "beta", "sigma"),
  labels = c("Loadings", "Means", "Residual SD"))]
pars.0.df.joined[, item := as.integer(gsub("[a-z]|\\[|\\]", "", variable))]
pars.0.df.joined

fig.capt <-
  "Indicator codes: 1 = Cooking fuel, 2 = Sanitation, 3 = Drinking water, 4 = Electricity, 5 = Housing, 6 = Assets"
ggplot(pars.0.df.joined, aes(item, m.1, fill = model, shape = model)) +
  geom_pointrange(aes(ymin = ll.1, ymax = ul.1), position = position_dodge(.5), alpha = .65) +
  theme_classic() + theme(legend.position = "top", strip.background = element_blank(),
                          axis.title.y = element_blank()) +
  scale_shape_manual(values = c(1, 4)) +
  facet_wrap(~ type, scales = "free_y") +
  labs(x = "Indicator number", subtitle = "Estimates with 95% interval",
       shape = "Method", fill = "Method", caption = fig.capt)
ggsave(paste0(print.images, "04_params_lav_bayes.pdf"), width = 6.5, height = 4)

rm(pars.0.df, pars.0.df.joined)
gc()

# Extract posterior of outcome variable
dim(yhat.0.f <- as.data.frame(fit.0, "yhat"))

yhat.0.df <- as.data.table(t(apply(yhat.0.f, 2, function (x) {
  c(mean = mean(x), var = var(x), quantile(x, c(.5, .25, .75, .025, .975)))
})))
yhat.0.df$item <- X$item
yhat.0.df$item_n <- X$item_n
yhat.0.df$id <- X$id
yhat.0.df$response <- X$response.l

setnames(yhat.0.df, c("50%", "25%", "75%", "2.5%", "97.5%"), c("median", "lli", "uli", "ll", "ul"))
yhat.0.df

yhat.0.df[, mean(response > ll & response < ul), item_n]
#            item_n        V1
# 1:   Cooking fuel 0.9800000
# 2:     Sanitation 0.9306931
# 3: Drinking water 0.9504950
# 4:    Electricity 1.0000000
# 5:        Housing 0.9500000
# 6:         Assets 0.9504950

yhat.0.df$theta <- apply(as.data.frame(fit.0, "theta_p"), 2, median)[yhat.0.df$id]

yhat.0.f.l <- melt(as.data.table(cbind(yhat.0.f, rep = 1:nrow(yhat.0.f))), id.vars = "rep")
yhat.0.f.l
yhat.0.f.l[, case := as.integer(gsub("yhat\\[|\\]", "", variable))]
yhat.0.f.l[, item := X$item[case]]
yhat.0.f.l[, id := X$id[case]]
yhat.0.f.l[, response := X$response.l[case]]
yhat.0.f.l[, resid := response - value]
yhat.0.f.l

yhat.0.f.w <- reshape(
  yhat.0.f.l[, .(id, item, rep, resid)], direction = "wide",
  v.names = c("resid"), idvar = c("id", "rep"), timevar = "item")
yhat.0.f.w <- yhat.0.f.w[order(rep, id)]
yhat.0.f.w

combn(6, 2)
yhat.0.f.cor.s <- cbind(
  yhat.0.f.w[, cor(resid.1, resid.2, use = "p"), by = rep], yhat.0.f.w[, cor(resid.1, resid.3, use = "p"), by = rep][, 2],
  yhat.0.f.w[, cor(resid.1, resid.4, use = "p"), by = rep][, 2], yhat.0.f.w[, cor(resid.1, resid.5, use = "p"), by = rep][, 2],
  yhat.0.f.w[, cor(resid.1, resid.6, use = "p"), by = rep][, 2], yhat.0.f.w[, cor(resid.2, resid.3, use = "p"), by = rep][, 2],
  yhat.0.f.w[, cor(resid.2, resid.4, use = "p"), by = rep][, 2], yhat.0.f.w[, cor(resid.2, resid.5, use = "p"), by = rep][, 2],
  yhat.0.f.w[, cor(resid.2, resid.6, use = "p"), by = rep][, 2], yhat.0.f.w[, cor(resid.3, resid.4, use = "p"), by = rep][, 2],
  yhat.0.f.w[, cor(resid.3, resid.5, use = "p"), by = rep][, 2], yhat.0.f.w[, cor(resid.3, resid.6, use = "p"), by = rep][, 2],
  yhat.0.f.w[, cor(resid.4, resid.5, use = "p"), by = rep][, 2], yhat.0.f.w[, cor(resid.4, resid.6, use = "p"), by = rep][, 2],
  yhat.0.f.w[, cor(resid.5, resid.6, use = "p"), by = rep][, 2]
)

# If wanting to test non-linear correlation, next 14 lines:
# yhat.0.f.w.na <- na.omit(yhat.0.f.w)
# library(energy)  # distance correlation measure
# yhat.0.f.w.na[rep < 4, dcor2d(resid.1, resid.2), rep]
# yhat.0.f.cor.s <- cbind(
#   yhat.0.f.w.na[, dcor2d(resid.1, resid.2), by = rep], yhat.0.f.w.na[, dcor2d(resid.1, resid.3), by = rep][, 2],
#   yhat.0.f.w.na[, dcor2d(resid.1, resid.4), by = rep][, 2], yhat.0.f.w.na[, dcor2d(resid.1, resid.5), by = rep][, 2],
#   yhat.0.f.w.na[, dcor2d(resid.1, resid.6), by = rep][, 2], yhat.0.f.w.na[, dcor2d(resid.2, resid.3), by = rep][, 2],
#   yhat.0.f.w.na[, dcor2d(resid.2, resid.4), by = rep][, 2], yhat.0.f.w.na[, dcor2d(resid.2, resid.5), by = rep][, 2],
#   yhat.0.f.w.na[, dcor2d(resid.2, resid.6), by = rep][, 2], yhat.0.f.w.na[, dcor2d(resid.3, resid.4), by = rep][, 2],
#   yhat.0.f.w.na[, dcor2d(resid.3, resid.5), by = rep][, 2], yhat.0.f.w.na[, dcor2d(resid.3, resid.6), by = rep][, 2],
#   yhat.0.f.w.na[, dcor2d(resid.4, resid.5), by = rep][, 2], yhat.0.f.w.na[, dcor2d(resid.4, resid.6), by = rep][, 2],
#   yhat.0.f.w.na[, dcor2d(resid.5, resid.6), by = rep][, 2]
# )

colnames(yhat.0.f.cor.s)[-1] <- paste0(
  "indicators: ", apply(combn(6, 2), 2, function (idxs) paste0(idxs, collapse = "-")))
yhat.0.f.cor.s

yhat.0.f.cor.s.l <- melt(yhat.0.f.cor.s, id.vars = "rep")
yhat.0.f.cor.s.l
yhat.0.f.cor.s.l <- yhat.0.f.cor.s.l[
  , .(median = median(value), iqr.lo = quantile(value, .25), iqr.hi = quantile(value, .75),
      lo = quantile(value, .025), hi = quantile(value, .975)), variable]
yhat.0.f.cor.s.l[, order := 1:.N]
yhat.0.f.cor.s.l
resid(lav.fit, type = "cor")$cov

ggplot(yhat.0.f.cor.s.l, aes(reorder(variable, -order), median)) +
  geom_linerange(aes(ymin = iqr.lo, ymax = iqr.hi), size = 2, alpha = .35) +
  geom_pointrange(aes(ymin = lo, ymax = hi)) +
  geom_hline(yintercept = 0, linetype = 2) + coord_flip() +
  theme(strip.background = element_blank(), panel.border = element_blank(),
        axis.ticks = element_blank()) +
  labs(y = "Correlation between residuals: median, IQR and 95% quantile intervals", x = "")
ggsave(paste0(print.images, "05_cor_resid_bayes.pdf"), width = 6.5, height = 2.5)

# Heteroskedastic pattern present (for cooking fuel and drinking water)
ggplot(yhat.0.df, aes(theta, response - median)) +
  geom_point(shape = 1, alpha = .25) +
  geom_abline(linetype = 2, intercept = 0, slope = 0) +
  facet_wrap(~ item_n) +
  geom_smooth(formula = y ~ s(x, bs = "ts"), method = "gam", col = cbbPalette[1],
              se = TRUE, size = .5) +
  geom_smooth(aes(y = response - ll), formula = y ~ s(x, bs = "ts"),
              method = "gam", col = cbbPalette[1], se = FALSE, size = .5) +
  geom_smooth(aes(y = response - ul), formula = y ~ s(x, bs = "ts"),
              method = "gam", col = cbbPalette[1], se = FALSE, size = .5) +
  theme_classic() + theme(strip.background = element_blank(), legend.position = "top") +
  labs(x = "Standardized country average (k_c^')", y = "residual with smoothed 95% interval", col = "Model")
ggsave(paste0(print.images, "06_resid_spread_bayes.pdf"), width = 6.5, height = 4.5)

# Compile Stan model (takes time)
mod.1 <- cmdstan_model(file.path(paste0(stan.scripts, "cong_joint.stan")))
fit.1 <- mod.1$sample(
  data = data,  # data list passed to Stan
  seed = 12345,  # Seed for reproducibility
  iter_warmup = 1e3,  # 1000 warmup samples per chain
  iter_sampling = 3e3,  # 3000 samples for inference per chain
  chains = 8,  # 8 chains
  parallel_chains = 8  # 8 cores on multicore systems
)
fit.1$cmdstan_diagnose()  # Run quick diagnostic checks
fit.1 <- read_stan_csv(fit.1$output_files())  # convert to Rstan

# Print major parameters
print(fit.1, c("Theta", "yhat", "lp__"), include = FALSE)
#                   mean se_mean   sd  2.5%   25%   50%   75% 97.5% n_eff Rhat
# alpha             0.55    0.01 0.79 -1.01  0.01  0.55  1.08  2.11  4408 1.00
# beta_dev[1]      -0.58    0.01 0.79 -2.13 -1.12 -0.57 -0.03  0.97  4591 1.00
# beta_dev[2]      -0.23    0.01 0.79 -1.78 -0.77 -0.22  0.31  1.30  4587 1.00
# beta_dev[3]       0.23    0.01 0.79 -1.33 -0.31  0.23  0.76  1.77  4567 1.00
# beta_dev[4]       0.30    0.01 0.79 -1.25 -0.24  0.30  0.84  1.84  4580 1.00
# beta_dev[5]      -0.16    0.01 0.79 -1.71 -0.70 -0.16  0.38  1.37  4577 1.00
# beta_dev[6]       0.75    0.01 0.79 -0.80  0.21  0.76  1.29  2.29  4590 1.00
# lambda[1]         1.28    0.00 0.11  1.08  1.20  1.27  1.35  1.50  2522 1.00
# lambda[2]         0.93    0.00 0.08  0.79  0.88  0.93  0.98  1.09  2463 1.00
# lambda[3]         0.62    0.00 0.06  0.51  0.58  0.62  0.66  0.74  3465 1.00
# lambda[4]         1.21    0.00 0.09  1.05  1.15  1.21  1.27  1.40  2119 1.00
# lambda[5]         0.87    0.00 0.08  0.72  0.81  0.86  0.92  1.02  3150 1.00
# lambda[6]         0.86    0.00 0.07  0.73  0.81  0.85  0.90  1.00  2544 1.00
# ln_alpha         -1.72    0.00 0.28 -2.27 -1.87 -1.72 -1.57 -1.17  3527 1.00
# ln_sigma2_dev[1]  0.61    0.01 0.33  0.00  0.40  0.60  0.81  1.29  3910 1.00
# ln_sigma2_dev[2] -0.07    0.01 0.33 -0.73 -0.27 -0.07  0.12  0.56  3835 1.00
# ln_sigma2_dev[3] -0.15    0.00 0.32 -0.78 -0.34 -0.15  0.04  0.47  4351 1.00
# ln_sigma2_dev[4] -0.58    0.00 0.36 -1.35 -0.80 -0.56 -0.35  0.05  5282 1.00
# ln_sigma2_dev[5]  0.39    0.00 0.31 -0.21  0.20  0.38  0.57  1.03  4007 1.00
# ln_sigma2_dev[6] -0.22    0.00 0.31 -0.84 -0.40 -0.22 -0.04  0.37  4175 1.00
# ls_scale          0.59    0.00 0.26  0.25  0.41  0.54  0.71  1.24  7255 1.00
# ll_scale          0.64    0.00 0.25  0.29  0.47  0.60  0.77  1.26 11353 1.00
# ln_lambda_pos[1]  0.68    0.00 0.24  0.25  0.52  0.67  0.83  1.19  6544 1.00
# ln_lambda_pos[2]  0.76    0.00 0.23  0.32  0.61  0.75  0.90  1.23  9174 1.00
# ln_lambda_oth[1] -0.47    0.01 0.32 -1.14 -0.67 -0.45 -0.25  0.09  3706 1.00
# ln_lambda_oth[2]  0.00    0.00 0.27 -0.53 -0.18  0.00  0.17  0.53  9304 1.00
# ln_lambda_oth[3] -0.33    0.00 0.27 -0.84 -0.50 -0.33 -0.16  0.21  7666 1.00
# ln_lambda_oth[4]  0.11    0.00 0.27 -0.39 -0.07  0.11  0.29  0.65  6320 1.00
# r_tr              0.82    0.00 0.07  0.67  0.78  0.83  0.87  0.95  1927 1.00
# ln_lambda[1]      0.68    0.00 0.24  0.25  0.52  0.67  0.83  1.19  6544 1.00
# ln_lambda[2]     -0.47    0.01 0.32 -1.14 -0.67 -0.45 -0.25  0.09  3706 1.00
# ln_lambda[3]      0.76    0.00 0.23  0.32  0.61  0.75  0.90  1.23  9174 1.00
# ln_lambda[4]      0.00    0.00 0.27 -0.53 -0.18  0.00  0.17  0.53  9304 1.00
# ln_lambda[5]     -0.33    0.00 0.27 -0.84 -0.50 -0.33 -0.16  0.21  7666 1.00
# ln_lambda[6]      0.11    0.00 0.27 -0.39 -0.07  0.11  0.29  0.65  6320 1.00
# beta[1]          -0.03    0.00 0.14 -0.31 -0.13 -0.03  0.07  0.25  1594 1.01
# beta[2]           0.31    0.00 0.10  0.11  0.25  0.32  0.38  0.51  1525 1.01
# beta[3]           0.77    0.00 0.08  0.63  0.72  0.77  0.83  0.92  1874 1.01
# beta[4]           0.84    0.00 0.13  0.59  0.76  0.84  0.93  1.08  1361 1.01
# beta[5]           0.38    0.00 0.10  0.18  0.32  0.38  0.45  0.58  1702 1.01
# beta[6]           1.30    0.00 0.09  1.11  1.24  1.30  1.36  1.48  1508 1.01
# sigma[1]          0.58    0.00 0.06  0.46  0.54  0.58  0.62  0.70  6946 1.00
# sigma[2]          0.41    0.00 0.04  0.33  0.38  0.41  0.44  0.49  7511 1.00
# sigma[3]          0.39    0.00 0.04  0.32  0.37  0.39  0.42  0.47  8855 1.00
# sigma[4]          0.32    0.00 0.04  0.24  0.29  0.32  0.34  0.40  9319 1.00
# sigma[5]          0.52    0.00 0.04  0.43  0.48  0.51  0.54  0.61 12316 1.00
# sigma[6]          0.38    0.00 0.03  0.32  0.36  0.38  0.40  0.45 21283 1.00
# r                 0.65    0.00 0.14  0.34  0.55  0.65  0.75  0.90  1927 1.00

# saveRDS(fit.1, "res/mod_2.rds")
fit.1 <- readRDS("res/mod_2.rds")
check_hmc_diagnostics(fit.1)
# 4 of 24000 iterations ended with a divergence (0.0166666666666667%).

# Traceplots:
mcmc_trace(fit.1, pars = c(
  paste0("beta[", 1:6, "]"), paste0("lambda[", 1:6, "]"),
  paste0("sigma[", 1:6, "]"), paste0("ln_lambda[", 1:6, "]")))
# Rank plots:
mcmc_rank_overlay(
  fit.1, pars = c(
    paste0("beta[", 1:6, "]"), paste0("lambda[", 1:6, "]"),
    paste0("sigma[", 1:6, "]"), paste0("ln_lambda[", 1:6, "]")),
  ref_line = TRUE)

{
  pars.1.df <- as.data.frame(summary(fit.1, c("ln_lambda"))$summary)
  pars.1.df$variable <- rownames(pars.1.df)
  pars.1.df
  pars.1.df <- as.data.table(pars.1.df[, c(4, 6, 8, 11)])
  pars.1.df
  setnames(pars.1.df, 1:3, c("ll", "median", "ul"))
  pars.1.df[, item := as.integer(gsub("[a-z]|\\[|\\]|_", "", variable))]
  pars.1.df
  pars.1.df$item_t <- levels(X$item_n)
  pars.1.df
}

plt.pars <- ggplot(pars.1.df, aes(reorder(item_t, -item), median)) + coord_flip() +
  geom_pointrange(aes(ymin = ll, ymax = ul), position = position_dodge(.5), alpha = .65, shape = 4) +
  theme_bw() + theme(legend.position = "top", strip.background = element_blank(),
                     axis.title.y = element_blank(), panel.border = element_blank(),
                     axis.ticks = element_blank()) +
  scale_shape_manual(values = c(1, 4)) +
  labs(x = "Indicator number", y = "Log-variance loadings with 95% quantile interval")
plt.pars

rm(pars.1.df)
gc()

dim(yhat.1.f <- as.data.frame(fit.1, "yhat"))

yhat.1.df <- as.data.table(t(apply(yhat.1.f, 2, quantile, c(.5, .25, .75, .025, .975))))
yhat.1.df$item <- X$item
yhat.1.df$item_n <- X$item_n
yhat.1.df$id <- X$id
yhat.1.df$response <- X$response.l

setnames(yhat.1.df, c("50%", "25%", "75%", "2.5%", "97.5%"), c("median", "lli", "uli", "ll", "ul"))
yhat.1.df

yhat.1.df[, mean(response > ll & response < ul), item_n]
#            item_n       V1
# 1:   Cooking fuel 0.980000
# 2:     Sanitation 0.990099
# 3: Drinking water 0.980198
# 4:    Electricity 1.000000
# 5:        Housing 0.950000
# 6:         Assets 0.960396

yhat.1.df$theta <- colMeans(as.data.frame(fit.1, "Theta"))[1:101][yhat.1.df$id]

# Heteroskedastic pattern present (for cooking fuel and drinking water)
plt.resid <- ggplot(yhat.1.df, aes(theta, response - median)) +
  geom_point(shape = 1, alpha = .25) +
  geom_abline(linetype = 2, intercept = 0, slope = 0) +
  facet_wrap(~ item_n, scales = "free_x") +
  geom_smooth(formula = y ~ s(x, bs = "ts"), method = "gam", col = cbbPalette[1],
              se = TRUE, size = .5) +
  geom_smooth(aes(y = response - ll), formula = y ~ s(x, bs = "ts"),
              method = "gam", col = cbbPalette[1], se = FALSE, size = .5) +
  geom_smooth(aes(y = response - ul), formula = y ~ s(x, bs = "ts"),
              method = "gam", col = cbbPalette[1], se = FALSE, size = .5) +
  theme_classic() + theme(strip.background = element_blank(), legend.position = "top") +
  labs(x = "Country averages", y = "residual with smoothed 95% interval", col = "Model")
(plt.pars + labs(tag = "A")) / (plt.resid + labs(tag = "B")) + plot_layout(heights = c(1.5, 4))
ggsave(paste0(print.images, "07_resid_bayes.pdf"), width = 6.5, height = 5.5)

yhat.1.df

ggplot(yhat.1.df, aes(theta, median)) +
  geom_point(aes(theta, response), shape = 1, alpha = .125, size = .5) +
  geom_ribbon(aes(ymin = ll, ymax = ul), alpha = .125, col = cbPalette[1]) +
  geom_line(alpha = .5, linetype = 2) +
  geom_smooth(aes(y = ll), formula = y ~ s(x, bs = "ts"), method = "gam",
              se = FALSE, size = .25, col = cbbPalette[1], linetype = 2) +
  geom_smooth(aes(y = ul), formula = y ~ s(x, bs = "ts"), method = "gam",
              se = FALSE, size = .25, col = cbbPalette[1], linetype = 2) +
  facet_wrap(item_n ~ ., scales = "fixed") +
  theme_classic() + theme(strip.background = element_blank(), legend.position = "top") +
  scale_color_manual(values = cbbPalette[-1]) +
  labs(x = "Standardized country averages", y = "95% quantile interval",
       col = "Model", linetype = "Model")
ggsave(paste0(print.images, "08_model_suggestions_1.pdf"), width = 6.5, height = 4)

ggplot(yhat.1.df, aes(theta, median)) +
  geom_point(aes(theta, response), shape = 1, alpha = .25, size = .5) +
  geom_line(alpha = .5) +
  geom_smooth(aes(y = ll), formula = y ~ s(x, bs = "ts"), method = "gam",
              se = FALSE, size = .5, col = cbbPalette[1]) +
  geom_smooth(aes(y = ul), formula = y ~ s(x, bs = "ts"), method = "gam",
              se = FALSE, size = .5, col = cbbPalette[1]) +
  facet_wrap(item_n ~ ., scales = "fixed", nrow = 1) +
  theme_classic() + theme(strip.background = element_blank(), legend.position = "top") +
  scale_color_manual(values = cbbPalette[-1]) +
  labs(x = "Standardized country averages", y = "Smoothed 95% quantile interval",
       col = "Model", linetype = "Model")
ggsave(paste0(print.images, "08_model_suggestions.pdf"), width = 6.5, height = 2.5)

Theta <- rbindlist(list(
  as.data.table(t(apply(as.data.frame(fit.0, "theta_p"), 2, function (x) {
    c(m = mean(x), s = sd(x), quantile(x, c(.05, .95)))
  }))),
  as.data.table(t(apply(as.data.frame(fit.1, "Theta"), 2, function (x) {
    c(m = mean(x), s = sd(x), quantile(x, c(.05, .95)))
  })))[1:101]
), idcol = "method")
Theta

Theta[, id := 1:.N, method]
Theta[, avg := mean(m), id]
Theta[, order := rank(avg)]
Theta[, method.t := c("Congeneric", "Location-Scale")[method]]
Theta

ggplot(Theta, aes(m, s)) +
  geom_smooth(se = FALSE, method = "gam", formula = y ~ s(x, bs = "ts"), col = cbbPalette[1], size = .5) +
  facet_wrap(~ method.t) + geom_point(shape = 1) +
  theme(strip.background = element_blank(), panel.border = element_blank(),
        panel.spacing = unit(.5, "cm"), axis.ticks = element_blank()) +
  labs(x = "Country averages (posterior mean)", y = "Country averages (posterior SD)")
ggsave(paste0(print.images, "09_posterior_mean.pdf"), width = 6.5, height = 4)

head(L.df <- as.matrix(as.data.frame(fit.0, "lambda")))
head(S.df <- as.matrix(as.data.frame(fit.0, "sigma")))

unique(sort(attr(lavPredict(lav.fit, se = "standard"), "se")[[1]])[1:10])
1 / sqrt(
  sum(lavInspect(lav.fit, "est")$lambda ^ 2 /
        diag(lavInspect(lav.fit, "est")$theta)))

hist(sem.s <- sapply(1:24e3, function (i) {
  1 / sqrt(sum(L.df[i, ] ^ 2 / S.df[i, ] ^ 2))
}), main = "", xlab = "")
print(describe(sem.s), digits = 4)
quantile(sem.s, c(.025, .5, .975))

(id.s <- rep(1:101, 6))
(item.s <- rep.int(1:6, rep(101, 6)))

head(t.df <- matrix(seq(-2, 2, length.out = 101)))
cov(t.df <- MASS::mvrnorm(101, rep(0, 2), matrix(c(1, 0, 0, 1), 2), empirical = TRUE))
describe(t.df)

head(theta1.df <- as.matrix(as.data.frame(fit.1, "Theta")[, 1:101]))
head(lm1.df <- as.matrix(as.data.frame(fit.1, "lambda")))
head(aln.df <- as.numeric(as.data.frame(fit.1, "ln_alpha")[, 1]))
head(b1ln.df <- as.matrix(as.data.frame(fit.1, "ln_sigma2_dev")))
head(lmln.df <- as.matrix(as.data.frame(fit.1, "ln_lambda")))
head(cor.df <- as.numeric(as.data.frame(fit.1, "r")[, 1]))

dim(T1.df <- sapply(1:24e3, function (i) {
  r <- cor.df[i]
  data.table(
    y = lm1.df[i, item.s] ^ 2 /
      exp(r * theta1.df[i, id.s] * lmln.df[i, item.s] + b1ln.df[i, item.s] + aln.df[i]),
    x = id.s)[, 1 / sqrt(sum(y)), x]$V1
}))
rm(lm1.df, aln.df, b1ln.df, lmln.df, cor.df)
gc()

T1.df.l <- as.data.table(melt(T1.df))
T1.df.l
T1.df.l <- T1.df.l[, .(median(value), quantile(value, .25), quantile(value, .75)), Var1]
T1.df.l[, id := colMeans(theta1.df)[Var1]]
T1.df.l
ggplot(T1.df.l, aes(id, V1)) +
  geom_line(alpha = .5, size = .5) +
  geom_line(aes(y = V2), alpha = .5, size = .5) +
  geom_line(aes(y = V3), alpha = .5, size = .5) +
  geom_hline(yintercept = quantile(sem.s, c(.5, .25, .75)), linetype = 2) +
  geom_rug(y = NA) +
  labs(y = "Standard error of measurement", x = "Location latent variable")

describe(Theta[method == 1, s])

lav.fit

X.mat <- as.matrix(dat[, 15:20])
X.mat
(p.means <- rowMeans(X.mat, na.rm = TRUE))
(missings <- which(is.na(X.mat), arr.ind = TRUE))
(X.mat[is.na(X.mat)] <- p.means[missings[, 1]])
X.mat <- t((t(X.mat) - colMeans(X.mat)) / apply(X.mat, 2, sd))
X.mat

(R.inv <- chol2inv(chol(cor(X.mat))))
(l.cor <- lavInspect(lav.fit, "est")$lambda / sqrt(
  lavInspect(lav.fit, "est")$lambda ^ 2 +
    diag(lavInspect(lav.fit, "est")$theta)
))
(R.inv %*% l.cor)
hist(X.mat %*% (R.inv %*% l.cor))
cor(X.mat %*% (R.inv %*% l.cor), lavPredict(lav.fit))
plot(X.mat %*% (R.inv %*% l.cor), lavPredict(lav.fit))

lambda0.df <- as.matrix(as.data.frame(fit.0, "lambda"))
sigma.df <- as.matrix(as.data.frame(fit.0, "sigma"))

dim(FS.dist <- sapply(1:24e3, function (i) {
  l <- lambda0.df[i, ]
  l <- l / sqrt(l ^ 2 + sigma.df[i, ] ^ 2)
  X.mat %*% (R.inv %*% l)
}))
describe(apply(FS.dist, 1, mean))
describe(apply(FS.dist, 1, sd))
plot(rowMeans(FS.dist), apply(FS.dist, 1, sd))

scatter.smooth(colMeans(as.data.frame(fit.0, "theta_p")), apply(as.data.frame(fit.0, "theta_p"), 2, sd))
scatter.smooth(apply(as.data.frame(fit.1, "Theta"), 2, sd)[1:101], apply(as.data.frame(fit.0, "theta_p"), 2, sd))
abline(a = 0, b = 1)
scatter.smooth(apply(as.data.frame(fit.1, "Theta"), 2, mean)[102:202],
               apply(as.data.frame(fit.1, "Theta"), 2, sd)[102:202])

# R Session Info (Contains exact package versions and system information) ----
sessionInfo()
# R version 4.0.3 (2020-10-10)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Pop!_OS 20.04 LTS
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/atlas/libblas.so.3.10.3
# LAPACK: /usr/lib/x86_64-linux-gnu/atlas/liblapack.so.3.10.3
# 
# locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
#  [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
#  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
# [10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#  [1] cmdstanr_0.2.2         rstan_2.21.3           StanHeaders_2.21.0-7   bayesplot_1.7.2       
#  [5] ggrepel_0.8.2          directlabels_2020.6.17 ggforce_0.3.2          ggplot2_3.3.3         
#  [9] latex2exp_0.4.0        patchwork_1.1.0        scales_1.1.1           semPlot_1.1.2         
# [13] lavaan_0.6-7           psych_2.0.9            data.table_1.14.0     
# 
# loaded via a namespace (and not attached):
#  [1] backports_1.2.1     Hmisc_4.4-2         BDgraph_2.63        plyr_1.8.6         
#  [5] igraph_1.2.6        splines_4.0.3       inline_0.3.17       digest_0.6.27      
#  [9] htmltools_0.5.1.1   matrixcalc_1.0-3    fansi_0.4.2         magrittr_2.0.1     
# [13] Rsolnp_1.16         checkmate_2.0.0     lisrelToR_0.1.4     cluster_2.1.0      
# [17] openxlsx_4.2.3      RcppParallel_5.0.2  matrixStats_0.58.0  prettyunits_1.1.1  
# [21] jpeg_0.1-8.1        sem_3.1-11          colorspace_2.0-1    xfun_0.19          
# [25] dplyr_1.0.5         callr_3.5.1         crayon_1.4.1        jsonlite_1.7.2     
# [29] lme4_1.1-26         regsem_1.6.2        survival_3.2-7      glue_1.4.2         
# [33] polyclip_1.10-0     gtable_0.3.0        mi_1.0              V8_3.4.0           
# [37] pkgbuild_1.2.0      abind_1.4-5         DBI_1.1.0           Rcpp_1.0.6         
# [41] xtable_1.8-4        htmlTable_2.1.0     tmvnsim_1.0-2       foreign_0.8-80     
# [45] Formula_1.2-4       stats4_4.0.3        truncnorm_1.0-8     htmlwidgets_1.5.3  
# [49] RColorBrewer_1.1-2  posterior_0.1.3     ellipsis_0.3.2      pkgconfig_2.0.3    
# [53] loo_2.4.1           XML_3.99-0.5        farver_2.1.0        nnet_7.3-14        
# [57] kutils_1.70         utf8_1.2.1          labeling_0.4.2      tidyselect_1.1.1   
# [61] rlang_0.4.11        reshape2_1.4.4      munsell_0.5.0       tools_4.0.3        
# [65] cli_2.5.0           generics_0.1.0      ggridges_0.5.2      fdrtool_1.2.15     
# [69] stringr_1.4.0       arm_1.11-2          processx_3.4.5      knitr_1.30         
# [73] zip_2.1.1           purrr_0.3.4         glasso_1.11         pbapply_1.4-3      
# [77] nlme_3.1-149        whisker_0.4         compiler_4.0.3      rstudioapi_0.13    
# [81] curl_4.3.1          png_0.1-7           huge_1.3.4.1        tibble_3.1.1       
# [85] statmod_1.4.35      tweenr_1.0.1        pbivnorm_0.6.0      stringi_1.5.3      
# [89] ps_1.5.0            qgraph_1.6.5        rockchalk_1.8.144   lattice_0.20-41    
# [93] Matrix_1.2-18       nloptr_1.2.2.2      vctrs_0.3.8         pillar_1.6.0       
# [97] lifecycle_1.0.0     OpenMx_2.18.1       corpcor_1.6.9       R6_2.5.0           
# [101] latticeExtra_0.6-29 gridExtra_2.3       codetools_0.2-16    boot_1.3-25        
# [105] MASS_7.3-53         gtools_3.8.2        assertthat_0.2.1    rjson_0.2.20       
# [109] withr_2.4.2         mnormt_2.0.2        mgcv_1.8-33         parallel_4.0.3     
# [113] quadprog_1.5-8      grid_4.0.3          rpart_4.1-15        coda_0.19-4        
# [117] minqa_1.2.4         carData_3.0-4       d3Network_0.5.2.1   base64enc_0.1-3    

