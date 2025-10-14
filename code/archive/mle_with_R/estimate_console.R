################################################################################
# Script: estimate.R 
# Author: Chris Campos 
#
# Estimates the demand model with parallel processing
# 
# Dependencies:
# (1) get_data.R
# (2) likelihood.R
# (3) objfn.R

################################################################################
# Set personal library path
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))

#install.packages(c("haven", "doParallel", "foreach", "numDeriv", "parallel"))

rm(list = ls())
gc()

win_os <- TRUE
dir <-  "/project/lausd/decentralized_choice"
setwd(dir)

if (!win_os){
  setwd("/Volumes/lausd/decentralized_choice")
  dir <- "/Volumes/lausd/decentralized_choice"
}

# Log file setup
log_file <- file.path(getwd(), "estimation_log_runs.txt")


# Load necessary libraries
library(parallel)    # For parallel processing
library(numDeriv)    # For numerical derivatives
library(haven)       # For reading Stata files
library(doParallel)  # For foreach parallel backend
library(foreach)     # For parallel loops

# Set up number of cores to use (adjust based on your machine)
num_cores <- detectCores() - 1


# Source helper files 
source(paste0(dir, "/code/get_data.R"))
source(paste0(dir, "/code/get_pi.R"))
source(paste0(dir, "/code/likelihood.R"))
source(paste0(dir, "/code/objfn.R"))

# Basic setup
J <- 43                     # Number of schools
R <- 5                    # Number of simulation draws for estimation
R_post <- 1000              # Number of draws for calculating posteriors
K <- 1                      # Number of mass points
lambda <- 0.05              # Smoothing parameter
homog <- 1                  # Treat charter schools as homogeneous (1) or heterogeneous (0)
parallel_flag <- 0          # Use parallel processing? (1 = yes, 0 = no)


# Maximization options (to be passed to optim)
options_fmin <- list(maxit = 500, reltol = 1e-6)

# Grade and subject for outcome model
grade <- 5
subject <- "math"



# Matrices of possible A and Z values
# Can only receive one offer 
Zmat <- diag(J)
Zmat <- rbind(Zmat, rep(0, J))
Amat <- Zmat

# Read raw data
rawdata <- read_dta(paste0(dir, "/data/structural_data_2004_2008_sample.dta"))
# Use the R version of get_data to extract variables:
data_list <- get_data(rawdata, J)


N <- data_list$N
A <- data_list$A
Achoice <- data_list$Achoice
A_dum <- data_list$Adummies  # Application choice dummies
Z <- data_list$Z
E <- data_list$E
Echoice <- data_list$Echoice
Y <- data_list$Y
X <- data_list$X
d <- data_list$d
yearapp <- data_list$yearapp
cohort <- data_list$cohort
lat <- data_list$lat
lon <- data_list$lon

# Reformat choice matrices:
Q_A <- nrow(Amat)
Q_Z <- nrow(Zmat)

# Useful for vectorized calculations later 
# Replicate Zmat into a 3D array: dimensions N x Q_Z x J
Z_big <- array(rep(Zmat, times = N), dim = c(N, Q_Z, J))
# Replicate Amat into a 3D array: dimensions N x Q_A x J
A_big <- array(rep(Amat, times = N), dim = c(N, Q_A, J))


# Get admission probabilities and the admission matrix:
pi_i <- as.matrix(read_dta(paste0(dir, "/data/student_p_i_2004_2008_sample.dta")))
pi_i <- pi_i[, 3:ncol(pi_i)]
storage.mode(pi_i) <- "numeric"

# Define explanatory variables
# Covariates (demeaned)
X <- scale(X, center = TRUE, scale = FALSE)
Q_X <- ncol(X)

# Distance variables:
X_2 <- cbind(1, X)  # Add an intercept
Q_D <- ncol(X_2)
# Create D as a 3D array with dimensions N x (Q_D+1) x J:
D_new <- array(0, dim = c(N, Q_D + 1, J))
for (j in 1:J) {
  # For the first Q_D columns, multiply each column of X_2 by d[,j]
  D_new[, 1:Q_D, j] <- X_2 * matrix(d[, j], nrow = N, ncol = Q_D, byrow = FALSE)
  # Last column is the squared distance:
  D_new[, Q_D + 1, j] <- d[, j]^2
}
D <- D_new
Q_D <- dim(D)[2]

# Application cost variables: use X_2 as W
W <- X_2
Q_W <- ncol(W)

# Initial values and simulation draws:
set.seed(20250627)  # For reproducibility
start_time <- Sys.time()

# We can parallelize the generation of simulation draws if needed
sims_theta <- array(rnorm(N * R), dim = c(N, R))
sims_c <- array(rnorm(N * R), dim = c(N, R))


# Set initial parameter vector
omega_choice_00 <- c(
  -7.11824914972133, -6.04424220264173, -6.40767817467762,
  -3.58927520448416,  3.90968340802561, -5.65532601463138,
  0.295578949038475,  1.94866293872075, -1.36364058560653,
  6.36804835063359,  -1.26396354596724,  4.94613747899377,
  3.99631310258744,   2.46752755803524,  -0.901547793912744,
  -3.86910452672889,   1.59013651121427,   0.770600385306727,
  -1.66666324726422,   0.876848201869338,   6.45980537291673,
  1.92234606355944, -11.5904929501861,    1.76605929904829,
  2.37299846624176,   7.22448170507261,   0.0692067143349845,
  -1.25334200133401,   2.14848173598695,   1.08698663579163,
  0.837491596522269,  0.280099708638706,  -0.517266949639412,
  -1.14558587693566,  -0.754412362439559,   5.48434513636614
) /10

omega_choice_0 <- c(
  0.303522446152143, -0.478481733839735, -0.731657152646574,
  -0.0283770849392046, 0.210841282923994, -0.368588497197534,
  0.151690821800793, 0.0281902957967152, -0.252742281887619,
  0.0891552008454534, -0.0985073334817339, -0.0444177338061709,
  0.396938685234604, 0.221123973272241, -0.0726396252110982,
  -0.311247862962208, 0.183379493717612, 0.00904610666482618,
  -0.0920140019692782, 0.238376665000752, 0.480466765756683,
  0.194620454259735, -0.97845691042508, 0.442053372587564,
  -0.22778537169495, 0.91813862355021, -0.0132724522108005,
  -0.193872322138262, 0.42345401001241, 0.10538027544863,
  0.0458236214298319, 0.118172453979119, -0.0609056150306577,
  -0.318385315179998, -0.199586684395466, 0.607853690376942
)/100

Nblocks <-1
cl <- NULL
# Full matrix to pass to the maximization routine:
data_max <- list(Nblocks = Nblocks, A = A, A_dum = A_dum, Achoice=Achoice, Z = Z, E=E, X = X, D = D, W = W,
                 A_big = A_big, Z_big = Z_big, pi_i = pi_i, lambda = lambda, homog = homog,
                 sims_theta = sims_theta, sims_c = sims_c, cl = cl, log_file = log_file)

# MAXIMIZE OBJECTIVE:

start_time <- Sys.time()




# Go in the right direction 
#optim_result1 <- optim(omega_choice_0,
#                       fn = function(x) { objfn(x, data_max)$Q },
#                       method = "Nelder-Mead", 
#                       control = list(maxit = 18000))

optim_result <- optim(omega_choice_0,
                      fn = function(x) { objfn(x, data_max, grad=FALSE)$Q },
                      gr = function(x) { objfn(x, data_max, grad=TRUE)$G },
                      method = "BFGS", control = options_fmin)

omega_choice_hat1 <- optim_result1$par

# Get standard errors:
cat("Calculating standard errors...\n")
start_time <- Sys.time()
res_obj <- objfn(omega_choice_hat1, data_nm, grad = TRUE)
Q_hat <- res_obj$Q
G <- res_obj$G
G_i <- res_obj$G_i

## This calculation can potentially be parallelized too (do later)
#fisher <- matrix(0, nrow = length(omega_choice_hat1), ncol = length(omega_choice_hat1))
#for (i in 1:N) {
#    fisher <- fisher + (1 / N) * (matrix(G_i[i, ], ncol = 1) %*% matrix(G_i[i, ], nrow = 1))
#}

# Save results based on specification
results_choice0 <- list(omega_choice_hat = omega_choice_hat1, Q_hat = Q_hat, G = G, G_i = G_i,
                        omega_choice_0 = omega_choice_0)

if (homog == 0 && K > 1) {
  saveRDS(results_choice0_bfgs, file = "results_choice_bfgs.rds")
} else if (homog == 1) {
  saveRDS(results_choice0_bfgs, file = "results_choice_homog_bfgs.rds")
} else if (homog == 0 && K == 1) {
  saveRDS(results_choice0_bfgs, file = "results_choice_K1_bfgs.rds")
}


log_message(paste("Standard error calculation complete. Elapsed time:", round(difftime(Sys.time(), start_time, units="mins"), 2), "minutes"))
log_message(paste("Q_hat:", Q_hat))
log_message(paste("Length of omega_choice_hat:", length(omega_choice_hat)))
log_message("Parameter estimates:")
log_message(capture.output(print(cbind(omega_choice_0, omega_choice_hat1, SE_choice))))


# Stop the cluster if it was created
if (!is.null(cl)) {
  stopCluster(cl)
}
