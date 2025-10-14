################################################################################
# Script: lr_tests.R
# Author: Chris Campos
#
# Makes a table with LR test results 
#
#
################################################################################
# Set personal library path
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
#install.packages(c("haven", "doParallel", "foreach", "numDeriv", "parallel"))

rm(list = ls())
gc()

dir <- "/Volumes/lausd/decentralized_choice"
# dir <- "Z:/decentralized_choice"
setwd(dir)



# Load necessary libraries
library(parallel)    # For parallel processing
library(numDeriv)    # For numerical derivatives
library(haven)       # For reading Stata files
library(doParallel)  # For foreach parallel backend
library(foreach)     # For parallel loops
library(matlib)

library(reticulate) # For loading numpy array

# Explicitly point to the Python binary from the module system
# Note that failing to specify python path can cause issue when reading NPZ file:
# reticulate won't convert the numpy array to matrix. 

#reticulate::install_python("3.12")

#use_python("C:/Users/admin/AppData/Local/r-reticulate/r-reticulate/pyenv/pyenv-win/versions/3.12.3/python.exe", required = TRUE) # For windows
# use_python("C:/Users/Administrator/AppData/Local/Programs/Python/Python312/python.exe", required = TRUE) # For windows
#use_python("/apps/python/3.10/3.10.9/bin/python3", required = TRUE) # For Mercury Cluster Python 3.10
 use_python("/usr/bin/python3", required = TRUE) # For Mac users?

# Print Python config to confirm
cat("Using Python from:", py_config()$python, "\n")
source(paste0(dir, "/code/helper/read_npz.R"))


################################################################################
############################ Standard Preparation ##############################
################################################################################

# Set up number of cores to use (adjust based on your machine)
num_cores <- parallelly::availableCores() - 1


# Source helper files 
source(paste0(dir, "/code/helper/get_data.R"))
source(paste0(dir, "/code/helper/get_pi.R"))
source(paste0(dir, "/code/helper/likelihood_mixture_types.R"))
source(paste0(dir, "/code/helper/names.R"))

# Basic setup
last_yr_2013  <- TRUE               # If true, uses 2004-2013 estimation results
J <- ifelse(last_yr_2013, 53, 40)   # Number of schools
R <- 300                            # Number of simulation draws for estimation
R_post <- 1000                      # Number of draws for calculating posteriors
lambda <- 0.05                      # Smoothing parameter
homog <- 0                          # Treat charter schools as homogeneous (1) or heterogeneous (0)
eta_out <- TRUE                     # Is cost = exp() + eta or exp(eta)?
parallel_flag <- 0                  # Partition data into multiple blocks for parallel processing? (1 = yes, 0 = no)
posterior    <- TRUE                # If false, only calculate the standard errors of the omega


last_yr <- ifelse(last_yr_2013, 2013, 2008)



# Maximization options (to be passed to optim)
options_fmin <- list(maxit = 100000, reltol = 1e-8)

# Grade and subject for outcome model
grade <- 5
subject <- "math"
eta <- "eta_out"
homog <- "hetero"

K <- 1 

npz_name <- paste0("in_mag_ind_", eta, "_results_choice_", homog, "_K", K,
                   ifelse(last_yr == 2013, "_2013", ""), # Are use using apps beyond year 2008?
                   "_bfgs.npz")

npz_path <- file.path(dir, "estimates", npz_name)
results1 <- read_npz(npz_path)


K <-2


npz_name <- paste0("in_mag_ind_", eta, "_results_choice_", homog, "_K", K,
                   ifelse(last_yr == 2013, "_2013", ""), # Are use using apps beyond year 2008?
                   "_bfgs.npz")

npz_path <- file.path(dir, "estimates", npz_name)
results2 <- read_npz(npz_path)


K <- 3 


npz_name <- paste0("in_mag_ind_", eta, "_results_choice_", homog, "_K", K,
                   ifelse(last_yr == 2013, "_2013", ""), # Are use using apps beyond year 2008?
                   "_bfgs.npz")

npz_path <- file.path(dir, "estimates", npz_name)
results3 <- read_npz(npz_path)


# Number of parameters
omega1 <- results1$omega_choice_hat$tolist()
omega2 <- results2$omega_choice_hat$tolist()
omega3 <- results3$omega_choice_hat$tolist()

# Number of obs 
Nobs <- results1$G_i$tolist()
Nobs <- nrow(do.call(rbind, lapply(Nobs, unlist)))

# negative log likelihood value 
ll1 <- as.numeric((-results1$Q_hat)$tolist())
ll2 <- as.numeric((-results2$Q_hat)$tolist())
ll3 <- as.numeric((-results3$Q_hat)$tolist())

bic1 <- as.numeric((-2 * ll1)$tolist()) + log(Nobs) * length(omega1)
bic2 <- as.numeric((-2 * ll2)$tolist()) + log(Nobs) * length(omega2)
bic3 <- as.numeric((-2 * ll3)$tolist()) + log(Nobs) * length(omega3)

# LR test is -2 ( ll_unrestricted -ll_restricted)
lr_2_1 <- as.numeric((-2 * (ll2 -ll1))$tolist() )
lr_3_2 <- as.numeric((-2 * (ll3 -ll2))$tolist() )


bic_2_1 <- bic1 - bic2
bic_3_2 <- bic2 - bic3

# Chi-squared p-value 
p_2_1 <- pchisq(lr_2_1, df = length(omega2) - length(omega1), lower.tail=FALSE)
p_3_2 <- pchisq(lr_3_2, df = length(omega3) - length(omega2), lower.tail=FALSE)




## ---------- Format helpers ----------
fmt_int   <- function(x) formatC(x, format = "d", big.mark = ",")
fmt_1dp   <- function(x) formatC(x, digits = 1, format = "f", big.mark = ",")
fmt_chi   <- function(x) formatC(x, digits = 1, format = "f", big.mark = ",")
fmt_p     <- function(p) ifelse(p < 0.001, ".000", sub("^0", ".", sprintf("%.3f", p)))

## ---------- Degrees of freedom ----------
k1 <- length(omega1); k2 <- length(omega2); k3 <- length(omega3)
df_21 <- k2 - k1
df_32 <- k3 - k2

## ---------- (Minor fix) p-values ----------
# upper tail directly:
p_2_1 <- pchisq(lr_2_1, df = df_21, lower.tail = FALSE)
p_3_2 <- pchisq(lr_3_2, df = df_32, lower.tail = FALSE)

## ---------- Labels & values ----------
npar    <- c(k1, k2, k3)
loglik  <- c(ll1, ll2, ll3)
bic     <- c(bic1, bic2, bic3)

# LR rows: only populate the relevant column (K=2 vs K=1 goes under K=2, etc.)
chi_row <- c("",
             paste0(fmt_chi(lr_2_1), " (", df_21, ")"),
             paste0(fmt_chi(lr_3_2), " (", df_32, ")"))
p_row   <- c("",
             fmt_p(p_2_1),
             fmt_p(p_3_2))

## ---------- Build LaTeX ----------
tab_lines <- c(
  "\\begin{tabular}{lccc}",
  "\\toprule",
  " & \\multicolumn{1}{c}{K=1} & \\multicolumn{1}{c}{K=2} & \\multicolumn{1}{c}{K=3}\\\\",
  "\\midrule",
  paste0("Number of parameters & ", fmt_int(npar[1]), " & ", fmt_int(npar[2]), " & ", fmt_int(npar[3]), " \\\\"),
  paste0("Log likelihood & ", fmt_1dp(loglik[1]), " & ", fmt_1dp(loglik[2]), " & ", fmt_1dp(loglik[3]), " \\\\"),
  paste0("BIC & ", fmt_1dp(bic[1]), " & ", fmt_1dp(bic[2]), " & ", fmt_1dp(bic[3]), " \\\\"),
  paste0("& & & \\\\"),
  "Likelihood ratio tests: & & & \\\\[-0.6em]",
  paste0("\\quad $\\chi^2$ statistic (df) &  & ", chi_row[2], " & ", chi_row[3], " \\\\"),
  paste0("\\quad $p$-value &  & ", p_row[2], " & ", p_row[3], " \\\\"),
  "\\bottomrule",
  "\\end{tabular}",
)

## ---------- Write to .tex ----------
out_dir  <- file.path(dir, "output/tables")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
out_tex  <- file.path(out_dir, "lr_table_K.tex")
writeLines(tab_lines, con = out_tex)



