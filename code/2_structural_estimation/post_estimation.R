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

dir <- "/project/lausd/decentralized_choice"
interdir <- "/code/chris/decentralized_choice/code" # this must be edited if running for main replication 
# dir <- "Z:/decentralized_choice"
setwd(dir)

# Log file setup
log_file <- paste0(dir, interdir, "/2_structural_estimation/logs", "/post_estimation_", format(Sys.time(), "%Y_%m_%d_%H_%M_%S"), ".txt")
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

library(reticulate) # For loading numpy array

# Explicitly point to the Python binary from the module system
# Note that failing to specify python path can cause issue when reading NPZ file:
# reticulate won't convert the numpy array to matrix. 

#reticulate::install_python("3.12")

#use_python("C:/Users/admin/AppData/Local/r-reticulate/r-reticulate/pyenv/pyenv-win/versions/3.12.3/python.exe", required = TRUE) # For windows
# use_python("C:/Users/Administrator/AppData/Local/Programs/Python/Python312/python.exe", required = TRUE) # For windows
use_python("/apps/python/3.10/3.10.9/bin/python3", required = TRUE) # For Mercury Cluster Python 3.10
# use_python("/usr/bin/python3", required = TRUE) # For Mac users?

# Print Python config to confirm
cat("Using Python from:", py_config()$python, "\n")
source(paste0(dir, interdir, "/helper/read_npz.R"))


################################################################################
############################ Standard Preparation ##############################
################################################################################

# Set up number of cores to use (adjust based on your machine)
num_cores <- parallelly::availableCores() - 1
log_message(paste("Using", num_cores, "cores for parallel processing"))


# Source helper files 
source(paste0(dir, interdir, "/helper/get_data.R"))
source(paste0(dir, interdir, "/helper/get_pi.R"))
source(paste0(dir, interdir, "/helper/likelihood_mixture_types.R"))
source(paste0(dir, interdir, "/helper/names.R"))

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

log_message(paste("Basic setup, lambda:", lambda))
log_message(paste("Basic setup, R:", R))

# Maximization options (to be passed to optim)
options_fmin <- list(maxit = 100000, reltol = 1e-8)

# Grade and subject for outcome model
grade <- 5
subject <- "math"

if (posterior){
  # Start parallel pool
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

  rawdata <- read_dta(paste0(dir, "/data/structural_data_2004_", last_yr, ".dta"))
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
  pi_i <- as.matrix(read_dta(paste0(dir, "/data/student_p_i_2004_", last_yr, ".dta")))
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
  
  #log_message(paste("Simulation draws completed in", round(difftime(Sys.time(), start_time, units="mins"), 2), "minutes"))
}

################################################################################
########### Read in different estimates and calculate posteriors ###############
################################################################################
post_estimation_calculations <- function(eta, homog, K, posterior, last_yr=2008){
  #results <- readRDS(paste0(dir, "/estimates/", 
  #                          eta, 
  #                          "_results_choice_", 
  #                          homog, 
  #                          "_K", K, "_bfgs.rds"))
  
  
  npz_name <- paste0("in_mag_ind_", eta, "_results_choice_", homog, "_K", K,
                     ifelse(last_yr == 2013, "_2013", ""), # Are use using apps beyond year 2008?
                     "_bfgs.npz")
  
  npz_path <- file.path(dir, "estimates", npz_name)
  results <- read_npz(npz_path)

  print(K)
  print(eta)
  print("Objfn Value:")
  print(results$Q_hat)
  
  ############################# Standard errror calculations ###################
  # Calculate variance-covariance of model estimates 
  #G_i <- results$G_i
  #G_i_list <- G_i$tolist()
  #G_i <- do.call(rbind, lapply(G_i_list, unlist))
  #class(G_i)
  G_i <- results$G_i
  #print("Gradients:")
  #print(results$G_i)
  # If G_i needs to be a matrix, convert it directly
  if(!is.matrix(G_i)) {
    G_i <- as.matrix(G_i)
  }
  omega_choice_hat1<- results$omega_choice_hat
  fisher <- matrix(0, nrow = length(omega_choice_hat1), ncol = length(omega_choice_hat1))
  N <- nrow(G_i)
  
  for (i in 1:N) {
    fisher <- fisher + (1 / N) * (matrix(G_i[i, ], ncol = 1) %*% matrix(G_i[i, ], nrow = 1))
  }
  inv_fisher <- inv(fisher)
  
  # Calculate the variance-covariance matrix of the transformed estimates,
  # taking the exponential of the random coefficient estimates 
  #omega <- results$omega_choice_hat 
  #omega <- omega$tolist()
  #omega <- do.call(rbind, lapply(omega, unlist))
  
  omega <- results$omega_choice_hat
  # Convert to vector if needed
  omega <- as.vector(omega)
  
  # Get idx of omega corresponding to log(sd(theta)) or log(sd(eta))
  omega_names <- results$omega_names
  #omega_names_list <- omega_names$tolist()
  #omega_names <- do.call(rbind, lapply(omega_names_list, unlist))
  log_idx <- which(grepl("log\\(", omega_names))
  print("Exponentiating following vars:")
  print(omega_names[log_idx])
  
  # Transform log(sd) to sd by exponentiating
  omega[log_idx] <- exp(omega[log_idx])
  a <- ifelse(last_yr==2013, 13,0)
  if(K==1){
    # Jacobian of the transformed matrix -- other elements are just themselves 
    J <- diag(length(omega))
    for (l_idx in log_idx){ 
      J[l_idx, l_idx] <- omega[l_idx] # Derivative of exp(x) is exp(x)
    }
    vcov_omega <- J %*% inv_fisher %*% t(J)  # Variance-covariance matrix of transofrmed omega
    se <- sqrt(diag((1/N)*vcov_omega) )
  }
  if(K==2){
    # Type probabilities 
    p <- (exp(omega[83+a]))/(1+ exp(omega[83+a]))  # Probability of type 1
    print("Alpha parsed:")
    print(omega[83+a])
    print("Probability parsed:")
    print(p)
    mu2 <- (-1)*(omega[84+a]*p)/(1-p)  # Mean of type 2
    J_p_alpha <- p*(1-p)  # Jacobian of the probability transformation
    J_mu2_mu1 <- (-1)*(p)/(1-p) # Derivative of mu2 with respect to alpha
    J_mu2_alpha <- (-1)*omega[84+a]*(p/(1-p)) # Derivative of mu2 with respect to mu1
    # Going to add p and mu2 to omega vector 
    omega <- c(omega, p, mu2)
    print(length(omega))
    
    # Jacobian of the transformed matrix -- other elements are just themselves
    J <- diag(length(omega)-2)
    J <- rbind(J, matrix(0, nrow=2, ncol=(length(omega)-2) ) )  # Add rows for p and mu2
    print(dim(J))
    
    # Derivative of exp(x) is exp(x)
    for (l_idx in log_idx){ 
      J[l_idx, l_idx] <- omega[l_idx] 
    }
    J[85+a,83+a] <- J_p_alpha  # Derivative of p with respect to alpha
    J[86+a,84+a] <- J_mu2_mu1  # Derivative of mu2 with respect to mu1
    J[86+a,83+a] <- J_mu2_alpha  # Derivative of mu2 with respect to alpha

    vcov_omega <- J %*% inv_fisher %*% t(J)  # Variance-covariance matrix of transofrmed omega
    se <- sqrt(diag((1/N)*vcov_omega) )
  }
  if(K==3){

    p1 <- (exp(omega[84+a]) )/( 1 + exp(omega[84+a]) + exp(omega[85+a]))  # Probability of type 1
    p2 <- (exp(omega[85+a]))/(1+ exp(omega[84+a]) + exp(omega[85+a])) # Probability of type 2
    p3 <- 1-p1-p2
    
    mu3 <- (-1)*(p1*omega[86+a] + p2*omega[87+a])/(p3)  # Mean of type 3
    J_p1_alpha1 <- p1*(1-p1)  # Jacobian of the probability transformation for type 1
    J_p2_alpha2 <- p2*(1-p2)  # Jacobian of the probability transformation for type 2
    J_p1_alpha2 <- -p1*p2
    J_p2_alpha1 <- -p1*p2
    J_mu3_mu1 <- -(exp(omega[84+a]))
    J_mu3_mu2 <- -(exp(omega[85+a]))
    J_mu3_alpha1 <- omega[86+a]*exp(omega[84+a])
    J_mu3_alpha2 <- omega[87+a]*exp(omega[85+a])
    
    # Going to add p1, p2, p3, and mu3 to omega vector
    omega <- c(omega, p1, p2, mu3)
    print(length(omega))
    
    # Jacobian of the transformed matrix -- other elements are just themselves
    J <- diag(length(omega)-3)
    J <- rbind(J, matrix(0, nrow=3, ncol=(length(omega)-3) ) )  # Add rows for p1, p2, mu3
    print(dim(J))
    
    # Derivative of exp(x) is exp(x)
    for (l_idx in log_idx){ 
      J[l_idx, l_idx] <- omega[l_idx] 
    }
    J[88+a,84+a] <- J_p1_alpha1  # Derivative of p1 with respect to alpha1
    J[88+a,85+a] <- J_p1_alpha2  # Derivative of p1 with respect to alpha2
    J[89+a,84+a] <- J_p2_alpha1  # Derivative of p2 with respect to alpha1
    J[89+a,85+a] <- J_p2_alpha2  # Derivative of p2 with respect to alpha2
    J[90+a,86+a] <- J_mu3_mu1  # Derivative of mu3 with respect to mu1
    J[90+a,87+a] <- J_mu3_mu2  # Derivative of mu3 with respect to mu2
    J[90+a,84+a] <- J_mu3_alpha1  # Derivative of mu3 with respect to alpha1
    J[90+a,85+a] <- J_mu3_alpha2  # Derivative of mu3 with respect to alpha2
    
    vcov_omega <- J %*% inv_fisher %*% t(J)  # Variance-covariance matrix of transofrmed omega
    se <- sqrt(diag((1/N)*vcov_omega) )
  }
  # Create output dataset in csv format 
  names <- param_names(K, last_yr)
  estimates <- cbind(names, omega, se )
  write.csv(estimates, file= paste0(dir, "/estimates/", 
                                      "in_mag_ind_",
                                      eta, 
                                      "_estimates_choice_", 
                                      homog, 
                                      "_K", K, "_", last_yr, ".csv"))
  
  ############################# Posterior calculations #########################
    
  if (posterior){
      # Full matrix to pass to the maximization routine:
      data_max <- list(K=K, Nblocks = Nblocks, A = A, A_dum = A_dum, Achoice=Achoice, Z = Z, E=E, X = X, D = D, W = W,
                       A_big = A_big, Z_big = Z_big, pi_i = pi_i, lambda = lambda, homog = homog,
                       choice_col = choice_col, E_col = E_col, eta_out = eta_out, num_cores = num_cores,
                       sims_theta = sims_theta, sims_c = sims_c, cl = cl, log_file = log_file)
      
    
      
      print("Estimating posteriors...")
      omega <- results$omega_choice_hat 
      omega <- as.vector(omega)
      #omega <- omega$tolist()
      #omega <- do.call(rbind, lapply(omega, unlist))
      
      post_values <- likelihood_K_mixture(omega=omega, data=data_max, posterior=TRUE)
      
      # Keep just the posteriors with identifiers 
      pi_i <- as.matrix(read_dta(paste0(dir, "/data/student_p_i_2004_", last_yr, ".dta")))
      posteriors <- cbind(pi_i[,1:2], post_values$posterior_mu, post_values$posterior_var)
      write.csv(posteriors, file = paste0(dir, "/estimates/", 
                                "in_mag_ind_",
                                eta, 
                                "_posteriors_", 
                                homog, 
                                "_K", K, "_", last_yr, ".csv"))
    }
}

if (last_yr == 2008){
  post_estimation_calculations(eta="eta_in", homog="hetero", K=1, posterior = posterior)
  post_estimation_calculations(eta="eta_out", homog="hetero", K=1, posterior = posterior)
  post_estimation_calculations(eta="eta_in", homog="hetero", K=2, posterior = posterior)
  post_estimation_calculations(eta="eta_out", homog="hetero", K=2, posterior = posterior)
  post_estimation_calculations(eta="eta_in", homog="hetero", K=3, posterior = posterior)
}

results <- post_estimation_calculations(eta="eta_out", 
                             homog="hetero", 
                             K=3, 
                             posterior = posterior, 
                             last_yr = last_yr)

