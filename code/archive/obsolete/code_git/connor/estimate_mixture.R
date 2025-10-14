################################################################################
# Script: estimate.R 
# Author: Chris Campos, Ryan Lee, Connor Fogal
#
# Estimates the demand model with parallel processing (parallel currently not working)
# Uses L-BFGS-B optimization method
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
setwd("/project/lausd/decentralized_choice")
dir <- "/project/lausd/decentralized_choice"

# Basic setup
J <- 40                     # Number of schools
R <- 300                    # Number of simulation draws for estimation
R_post <- 1000              # Number of draws for calculating posteriors
K <- 2                      # Number of mass points
lambda <- 0.05              # Smoothing parameter
homog <- 0                  # Treat charter schools as homogeneous (1) or heterogeneous (0)
eta_out <- FALSE             # Is cost = exp() + eta or  exp(eta)?
parallel_flag <- 0          # Use parallel processing? (1 = yes, 0 = no) (Partition data into multiple blocks)
memoize       <- TRUE       # For calling objfn only once / iteration
                            # Need to study what's the risk of memoization. Reduction in time isn't large

# Log file setup
f_name <- paste("/eta", ifelse(eta_out, "out", "in"), "K", K, "start", sep = "_")
log_file <- paste0(getwd(), f_name, format(Sys.time(), "%Y_%m_%d_%H_%M_%S"), ".txt")
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

if (memoize){
  log_message("Using memoizaiton")
}

# Load necessary libraries
library(parallel)    # For parallel processing
library(numDeriv)    # For numerical derivatives
library(haven)       # For reading Stata files
library(doParallel)  # For foreach parallel backend
library(foreach)     # For parallel loops
library(memoise)

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

# Read raw data
log_message("Reading raw data...")
rawdata <- read_dta(paste0(dir, "/data/structural_data_2004_2008_sample.dta"))
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
set.seed(12345)  # For reproducibility
log_message("Generating simulation draws...")
start_time <- Sys.time()

# We can parallelize the generation of simulation draws if needed

sims_theta <- array(rnorm(N * R), dim = c(N, R))
sims_c <- array(rnorm(N * J* R), dim = c(N, J, R))

log_message(paste("Simulation draws completed in", round(difftime(Sys.time(), start_time, units="mins"), 2), "minutes"))
  
# Set initial parameter vector
set.seed(123456)
#omega_choice_0 <- rnorm(n=42, mean=0,sd=0.02)
#omega_choice_0[1] <- -1.4
#omega_choice_0[14] <- -0.4
#omega_choice_0[28] <- 0.2
#omega_choice_0[29] <- -0.1
#omega_choice_0[28] <- 0.1

omega_choice_0 <- c(
  -1.4567156534,  0.0094458423, -0.4045561186,  0.0598194427,  0.1262253800,
  0.0444772630, -0.0246325915,  0.0127847843,  0.0204510921,  0.0572385747,
  0.0399978301,  0.1046140873, -0.0484484979, -0.2616026753, -0.0067061133,
  0.0436389238,  0.0065852441,  0.0096586141, -0.0039726424, -0.0031864606,
  0.0121617194,  0.0033748398, -0.0034557826, -0.0148961058,  0.1164555316,
  0.0287177805,  0.0006761936, -0.0933390767, -1.2548835300, -0.1837808304,
  0.0150895013, -0.1734069139,  0.4206879331,  0.1310351966,  0.0523076142,
  0.1489851766,  0.0369591813, -0.2242043754, -0.1230576988,  0.0306908241,
  -1.2300950404, -0.3488746794
)

if (homog == 0){
  
  if (eta_out){
    
    omega_choice_0 <- c(
      -1.3902371703, -1.8526023375,  0.3009050067,  0.1998555338, -2.9830672150,
      -1.4328040076, -2.5842449671, -2.9175639202, -3.3778959913, -1.0457194870,
      -2.2269435274, -1.6815559729, -1.7896967035,  0.3932910146, -3.4813256170,
      -1.8548114638, -2.5444085078, -0.7304389022, -0.3855495949, -3.0635995822,
      -1.9992970561, -1.8647098574, -0.4580687829, -2.5111326784, -1.3014435751,
      -2.7304332905, -1.2855214072, -3.4199850531, -3.3506015724,  0.7493261333,
      -2.1805915698, -2.3581213990, -1.4151000380, -0.8281221977,  1.2769356630,
      -2.2037145933, -2.7442227473, -1.4026167491,  1.5170770975,  7.6319743969,
      -0.0192539655, -0.0957056704, -0.2241712754, -1.1580334301,  0.1020350334,
      -0.3270059392, -0.0661048711,  0.2483186936,  0.1990287496, -0.0660202129,
      1.5564785118,  0.1608986388,
      -0.2832288332, -0.0056573755,  0.0084257195, -0.0129952532, -0.1695916140,
      -0.0215097174,  0.0043080743,  0.0087984516,  0.0012777543, -0.0082398183,
      -0.0008711711,  0.0689624232,  0.0245605557,
      0.0027716644, -0.2784696062, -0.0034393270, -0.0574381677,  0.0644226877,
      -0.1147370817, -0.0012846066, -0.0014528080, -0.0173120877, -0.0363893057,
      -0.0253497484,  0.0074928664, -0.0892310992, -0.0355494132,
      -1.5609554828, -0.2373604967, -0.1076556340,
      2.5353127870, -0.2317632725
    )
    
    
  } else{
    # Vector where most recent run crashed at 
    omega_choice_0 <- c(
      -0.413623808, -2.580312098,  4.319178899,  2.508565641, -4.817763414,
      -0.716621444, -3.480940590, -4.117435268, -5.398003961,  0.021385933,
      -2.674465131, -0.969820811, -0.216500942,  3.699898157, -6.034771722,
      -1.356530864, -3.411124616,  0.508582554,  2.231623942, -4.331834160,
      -1.703126039, -1.929706788,  1.333753763, -2.656508963,  0.460670680,
      -3.477386189,  0.934638905, -5.639295327, -4.758858777,  4.635361207,
      -1.820286695, -3.268209320, -0.792099335,  1.195799868,  5.180838665,
      -1.539783918, -4.178581642, -0.163993437,  5.482039769, 12.936861816,
      -0.154073914,  0.586549266,  0.278604048, -0.800198342,  0.506130207,
      -0.810743705, -0.014257831,  0.477572693,  0.521899408, -0.107502877,
      3.632301639, -0.223510608,
      -0.817777561, -0.026221998, -0.009170093, -0.070483907, -0.246081600,
      -0.010398912,  0.033353882, -0.004932554, -0.002599460, -0.005359491,
      -0.022438172,  0.063568884,  0.052739429,
      0.014604937,  2.651679538, -0.093032393,  0.057889442,  0.508409183,
      0.047089518,  0.111378576,  0.088785725, -0.081880660, -0.161952718,
      -0.087403626,  0.008733317, -0.516750831, -0.297260301,
      0.480225128,  0.501874224,  0.586712173,
      2.337216832, -0.670004060
    )
    
  }

  
  
  
}

print("starting values:")
print(omega_choice_0)


# Full matrix to pass to the maximization routine:
data_max <- list(K=K, Nblocks = Nblocks, A = A, A_dum = A_dum, Achoice=Achoice, Z = Z, E=E, X = X, D = D, W = W,
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


# MAXIMIZE OBJECTIVE:
#log_message("Starting Nelder-Mead optimization...")
# start_time <- Sys.time()
# Go in the right direction 
#optim_result1 <- optim(omega_choice_0,
#                       fn = function(x) { objfn_no_partition(x, data_max)$Q },
#                       method = "Nelder-Mead", 
#                       control = list(maxit = 3))
#omega_choice_1 <- optim_result1$par
#start_time <- Sys.time()
#objfn_no_partition(omega_choice_0, data_max, grad=TRUE)
#print(Sys.time()-start_time) 



print("starting values:")
print(omega_choice_0)

#log_message("Starting NM optimization...")
#optim_result_nm <- optim(omega_choice_0,
#                         fn = function(x) {objfn_no_partition(x, data_max, grad=FALSE)$Q },
#                         control = options_fmin)
#omega_choice_nm <- optim_result_nm$par


log_message("Starting BFGS optimization...")

if (!memoize){
  # Maximization options (to be passed to optim)
  optim_result <- optim(omega_choice_0,
                        fn = function(x) { objfn_no_partition(x, data_max, grad=FALSE)$Q },
                        gr = function(x) { objfn_no_partition(x, data_max, grad=TRUE)$G },
                        method = "L-BFGS-B", 
                        lower  = -Inf,
                        upper  =  Inf,
                        control = options_fmin)
} else{
  objfn_wrap <- function(par){
    objfn_no_partition(omega=par, data=data_max, grad=TRUE) 
  }
  
  mlmem <- memoise(objfn_wrap)
  
  fn_mem <- function(par){
    mlmem(par)$Q
  }
  gr_mem <- function(par){
    mlmem(par)$G
  }

  optim_result <- optim(omega_choice_0,
                        fn = fn_mem,
                        gr = gr_mem,
                        method = "L-BFGS-B", 
                        lower  = -Inf,
                        upper  =  Inf,
                        control = options_fmin)
}

omega_choice_hat1 <- optim_result$par



# Get standard errors:
start_time <- Sys.time()
res_obj <- objfn_no_partition(omega_choice_hat1, data_max, grad = TRUE)
Q_hat <- res_obj$Q
G <- res_obj$G
G_i <- res_obj$G_i

# Print out the optimization results (e.g. convergence indicator)
print(optim_result)

# Save results based on specification
results_choice0_bfgs <- list(omega_choice_hat = omega_choice_hat1, 
                             Q_hat = Q_hat, 
                             G = G, G_i = G_i,
                             omega_choice_0 = omega_choice_0, 
                             optim_result=optim_result, 
                              data_max = data_max, log_file = log_file
                             )

if (homog == 0) {
  saveRDS(results_choice0_bfgs, file = paste0("eta_", ifelse(eta_out, "out", "in"), "_results_choice_K", K, "_bfgs.rds"))
} else if (homog == 1) {
  saveRDS(results_choice0_bfgs, file = "results_choice_homog_bfgs.rds")
} 

log_message(paste("Standard error calculation complete. Elapsed time:", round(difftime(Sys.time(), start_time, units="mins"), 2), "minutes"))
log_message(paste("Q_hat:", Q_hat))
log_message(paste("Length of omega_choice_hat:", length(omega_choice_hat1)))
log_message("Parameter estimates:")
log_message(capture.output(print(cbind(omega_choice_0, omega_choice_hat1))))

# Stop the cluster if it was created
if (!is.null(cl)) {
  stopCluster(cl)
}
