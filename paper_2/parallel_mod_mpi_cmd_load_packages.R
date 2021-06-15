# Package notes ----

# Installation steps for CommandStan:
# https://github.com/stan-dev/cmdstan/wiki/Getting-Started-with-CmdStan

# In addition to the packages list below, users should install:
# - posterior: install.packages("posterior", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

# Load packages
library(data.table)  # Fast alternative to data frames
library(psych)  # Dataset descriptives
library(lavaan)  # Standard SEM
library(semPlot)  # plot lavaan fit
library(scales)  # Formatting outputs
library(patchwork)  # Adding ggplots together
library(latex2exp)  # Creating equations using latex
library(ggplot2)  # Plotting
library(ggforce)  # extend ggplot
library(directlabels)  # Direct labeling of ggplots
library(ggrepel)  # Labeling of lines on ggplot
theme_set(theme_bw())  # ggplot label
library(ggforce)  # Enhance ggplot2 functionalities

# Load Stan
library(bayesplot)  # Plotting aspects of Stan fit
color_scheme_set("viridis")  # Color scheme for plots
library(rstan)  # Load Rstan
library(cmdstanr)  # Command Stan
set_cmdstan_path("~/cmdstan/")  # Command Stan path
