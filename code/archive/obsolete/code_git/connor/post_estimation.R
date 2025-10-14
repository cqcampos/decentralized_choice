################################################################################
# Script: post_estimation.R
# Author: Chris Campos
#
# Takes demand estimates as given and does the following:
# 1. Calculates standard errors, using delta method where appropriate 
# 2. Makes a table with estimates and standard errors 
# 3. Calculates posterior means of theta for each observation 
# 4. Calculates the preference index, observable, unobservable, and combined
# 5. Outputs a dataset with posteriors and the preference indices with studentpseudo links
#
#
################################################################################
# Set personal library path
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
#install.packages(c("haven", "doParallel", "foreach", "numDeriv", "parallel"))

rm(list = ls())
gc()
setwd("/project/lausd/decentralized_choice")
dir <- "/project/lausd/decentralized_choice"

# Log file setup
log_file <- paste0(getwd(), "/post_estimation_", format(Sys.time(), "%Y_%m_%d_%H_%M_%S"), ".txt")
# Initialize log file with timestamp
write(paste("=== Estimation started at", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "==="), 
      file = log_file)

# Custom logging function that writes to both console and log file
log_message <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  formatted_msg <- paste0("[", timestamp, "] ", msg)
  # Print to console
  cat(formatted_msg, "\n")
  # Append to log file
  write(formatted_msg, file = log_file, append = TRUE)
}

log_message("Starting estimation process...")

# Load necessary libraries
library(parallel)    # For parallel processing
library(numDeriv)    # For numerical derivatives
library(haven)       # For reading Stata files
library(doParallel)  # For foreach parallel backend
library(foreach)     # For parallel loops
library(matlib)


################################################################################
############################ Standard Preparation ##############################
################################################################################

# Set up number of cores to use (adjust based on your machine)
num_cores <- parallelly::availableCores() - 1
log_message(paste("Using", num_cores, "cores for parallel processing"))


# Source helper files 
source(paste0(dir, "/code/get_data.R"))
source(paste0(dir, "/code/get_pi.R"))
source(paste0(dir, "/code/likelihood.R"))
source(paste0(dir, "/code/objfn.R"))
source(paste0(dir, "/code/objfn_no_partition.R"))
source(paste0(dir, "/code/likelihood_mixture_types.R"))
source(paste0(dir, "/code/names.R"))

# Basic setup
J <- 40                     # Number of schools
R <- 300                   # Number of simulation draws for estimation
R_post <- 1000              # Number of draws for calculating posteriors
K <- 1                      # Number of mass points
lambda <- 0.05              # Smoothing parameter
homog <- 0                  # Treat charter schools as homogeneous (1) or heterogeneous (0)
eta_out <- FALSE            # Is cost = exp() + eta or  exp(eta)?
parallel_flag <- 0          # Use parallel processing? (1 = yes, 0 = no) (Partition data into multiple blocks)

log_message(paste("Basic setup, lambda:", lambda))
log_message(paste("Basic setup, R:", R))

# Maximization options (to be passed to optim)
options_fmin <- list(maxit = 100000, reltol = 1e-8)

# Grade and subject for outcome model
grade <- 5
subject <- "math"

# Start parallel pool
# Ryan edit: I edited objfn, so that boolean `use_parallel` is = 1 if Nblocks > 1
#             Without this change, fact that !is.null(cl) will make `use_parallel` = 1, which we don't want
# cl <- makeCluster(num_cores)
cl <- NULL
Nblocks <- 1

if (parallel_flag == 1) {
  Nblocks <- num_cores 
  cl <- makeCluster(Nblocks)
  registerDoParallel(cl)
}

# Matrices of possible A and Z values
# Can only receive one offer 
Zmat <- diag(J)
Zmat <- rbind(Zmat, rep(0, J))
Amat <- Zmat


################################################################################
################### Read in full sample and data prep  #########################
################################################################################
# Read raw data
log_message("Reading raw data...")
rawdata <- read_dta(paste0(dir, "/data/structural_data_2004_2008.dta"))
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
choice_col <- data_list$choice_col
E_col <- data_list$E_col

log_message(paste("Data loaded. Number of observations:", N))

# Reformat choice matrices:
Q_A <- nrow(Amat)
Q_Z <- nrow(Zmat)

# Useful for vectorized calculations later 
# Replicate Zmat into a 3D array: dimensions N x Q_Z x J
Z_big <- array(rep(Zmat, times = N), dim = c(N, Q_Z, J))
# Replicate Amat into a 3D array: dimensions N x Q_A x J
A_big <- array(rep(Amat, times = N), dim = c(N, Q_A, J))


# Get admission probabilities and the admission matrix:
pi_i <- as.matrix(read_dta(paste0(dir, "/data/student_p_i_2004_2008.dta")))
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
set.seed(12345)  # For reproducibility
log_message("Generating simulation draws...")
start_time <- Sys.time()

sims_theta <- array(rnorm(N * R), dim = c(N, R))
sims_c <- array(rnorm(N * J* R), dim = c(N, J, R))

log_message(paste("Simulation draws completed in", round(difftime(Sys.time(), start_time, units="mins"), 2), "minutes"))

################################################################################
########### Read in different estimates and calculate posteriors ###############
################################################################################
post_estimation_calculations <- function(eta, homog, K){
  results <- readRDS(paste0(dir, "/estimates/", 
                            eta, 
                            "_results_choice_", 
                            homog, 
                            "_K", K, "_bfgs.rds"))
  ############################# Standard errror calculations ###################
  # Calculate variance-covariance of model estimates 
  G_i <- results$G_i
  omega_choice_hat1<- results$omega_choice_hat
  fisher <- matrix(0, nrow = length(omega_choice_hat1), ncol = length(omega_choice_hat1))
  N <- nrow(G_i)
  for (i in 1:N) {
    fisher <- fisher + (1 / N) * (matrix(G_i[i, ], ncol = 1) %*% matrix(G_i[i, ], nrow = 1))
  }
  inv_fisher <- inv(fisher)
  
  # Calculate the variance-covariance matrix of the transformed estimates,
  # taking the exponential of the random coefficient estimates 
  omega <- results$omega_choice_hat 
  omega[67:68] <- exp(omega[67:68])  # Exponentiate the random coefficients
  # Jacobian of the transformed matrix -- other elements are just themselves 
  J <- diag(length(omega))
  J[67, 67] <- diag(exp(omega[67]))  # Derivative of exp(x) is exp(x)
  J[68, 68] <- diag(exp(omega[68]))  # Derivative of exp(x) is exp(x)
  se <- sqrt(diag((1/N)*vcov_omega) )
  # Create output dataset in csv format 
  names <- param_names(K)
  estimates <- cbind(names, omega, se )
  write.csv(estimates, file= paste0(dir, "/estimates/", 
                                    eta, 
                                    "_estimates_choice_", 
                                    homog, 
                                    "_K", K, ".csv"))
  
  ############################# Posterior calculations #########################
  
  # Full matrix to pass to the maximization routine:
  data_max <- list(Nblocks = Nblocks, A = A, A_dum = A_dum, Achoice=Achoice, Z = Z, E=E, X = X, D = D, W = W,
                   A_big = A_big, Z_big = Z_big, pi_i = pi_i, lambda = lambda, homog = homog,
                   choice_col = choice_col, E_col = E_col, eta_out = eta_out, num_cores = num_cores,
                   sims_theta = sims_theta, sims_c = sims_c, cl = cl, log_file = log_file)
  
  # Remove redundant objects to save memory
  rm(data_list)
  rm(X)
  rm(Z)
  rm(E)
  rm(W)
  rm(A_dum)
  rm(A)
  rm(sims_theta)
  rm(sims_c)
  gc()
  
  print("Estimating posteriors")
  
  post_values <- likelihood_K_mixture(omega=results$omega_choice_hat, data=data_max, posterior=TRUE)
  
  # Keep just the posteriors with identifiers 
  pi_i <- as.matrix(read_dta(paste0(dir, "/data/student_p_i_2004_2008.dta")))
  posteriors <- cbind(pi_i[,1], post_values$posterior_mu)
  write.csv(posteriors, file = paste0(dir, "/estimates/", 
                                      eta, 
                                      "_posteriors_", 
                                      homog, 
                                      "_K", K, ".csv"))
}

post_estimation_calculations(eta="eta_in", homog="hetero", K=1)
post_estimation_calculations(eta="eta_out", homog="hetero", K=1)