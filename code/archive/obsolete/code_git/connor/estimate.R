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

# Log file setup
log_file <- paste0(getwd(), "/estimation_log_runs_lfbgs_", format(Sys.time(), "%Y_%m_%d_%H_%M_%S"), ".txt")
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

# Set up number of cores to use (adjust based on your machine)
num_cores <- parallelly::availableCores() - 1
log_message(paste("Using", num_cores, "cores for parallel processing"))


# Source helper files 
source(paste0(dir, "/code/get_data.R"))
source(paste0(dir, "/code/get_pi.R"))
source(paste0(dir, "/code/likelihood.R"))
source(paste0(dir, "/code/objfn.R"))
source(paste0(dir, "/code/objfn_no_partition.R"))

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
#cl <- makeCluster(num_cores)
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
  
  mean_utilities <- c(
    -1.5306461250, -1.4536306362, -1.4513813413, -1.4501382402,
    -1.5002890611, -1.4421260907, -1.4592557851, -1.4856363618, -1.4950136388,
    -1.4220915119, -1.4628642512, -1.4532892842, -1.4718431553, -1.4021597731,
    -1.4779642214, -1.4531418035, -1.4641341756, -1.4419304140, -1.4386325923,
    -1.4732622877, -1.4397752845, -1.4332621146, -1.4273440086, -1.4489139084,
    -1.4465724978, -1.4821651867, -1.4567884156, -1.4833180479, -1.4884737881,
    -1.3654952571, -1.4436104888, -1.4786891243, -1.4469822973, -1.4333286893,
    -1.4304179314, -1.4760034539, -1.4622433229, -1.4418455864, -1.3883943050,
    -1.4254272235
  )
  
  omega_choice_0 <- c(
    c(mean_utilities, omega_choice_0[-1])
  )
  
  
  omega_choice_0 <- c(
    -1.5273070486, -1.4505439182, -1.4478058585, -1.4485347921, -1.5251124148,
    -1.4297778675, -1.4625022770, -1.5038904933, -1.5161773503, -1.4049879499,
    -1.4669138843, -1.4466513744, -1.4868724399, -1.3801718574, -1.4973830874,
    -1.4495378590, -1.4647392202, -1.4355073088, -1.4329171392, -1.4842348071,
    -1.4302073427, -1.4206378548, -1.4148053360, -1.4512447987, -1.4411492403,
    -1.4966723488, -1.4573415382, -1.5008486150, -1.5053960716, -1.3218722748,
    -1.4489900965, -1.4839062698, -1.4421928343, -1.4226483776, -1.4212409552,
    -1.4990745354, -1.4676032302, -1.4338357875, -1.3674641352, -1.4248772634,
    0.0093429083, -0.4123784207,  0.0541479148,  0.1225125266,  0.0386712709,
    -0.0338734225,  0.0134439946,  0.0621339026,  0.1015327495,  0.1458073280,
    0.0933844803, -0.0313330093,
    -0.3575403233, -0.0033952946,  0.0485020910,  0.0109991609, -0.0050052987,
    -0.0131576610, -0.0004627495,  0.0251614963,  0.0074106864, -0.0063632248,
    -0.0258110419,  0.0772294740,  0.0486893627,
    0.0029949208, -0.1086637997, -1.1938037002,
    -0.1473563997,  0.0040118730, -0.1913924901,  0.3456318461,  0.2366780258,
    -0.0008626554,  0.1310154661,  0.0572712962, -0.1753681572, -0.0983553656,
    0.0336753695, -1.1030335759, -0.1962768795
  )
  
  # Vector where most recent run crashed at 
  omega_choice_0 <- c(
    -1.443764222, -1.736568081,  0.429236463, -1.023745931, -2.174548727,
    -1.041304409, -1.834764331, -2.138159699, -2.355993973, -0.628297020,
    -1.983441614, -1.218230504, -1.700241086,  0.460784801, -2.316723068,
    -2.510934918, -1.667738996, -0.869248333, -0.618294190, -2.157978229,
    -1.372681046, -1.256844249, -0.636308485, -1.734697852, -1.257531486,
    -1.959322014, -1.189485071, -2.306680023, -2.729375591,  0.688615641,
    -1.900825960, -2.041522274, -1.201060784, -0.892972511,  0.823543569,
    -1.964987274, -1.956701128, -1.009882851,  1.578319385, -10.617364775,
    -0.087884213,  1.032798370,  0.785773915,  0.104564089,  0.328069601,
    -0.067432694, -0.134315052,  0.218182862,  0.261970913, -0.017758676,
    0.661437984,  0.606130571,
    -0.309198844, -0.008358335, -0.059917857, -0.049511402, -0.195251183,
    -0.003604148,  0.005287851,  0.027569121,  0.019936418, -0.019221638,
    -0.006308082,  0.136915884,  0.010118429,
    0.001617983, -0.240497101,  0.105596238,
    0.268214192, -0.070057654,  0.169699020,  0.476557261, -0.074515275,
    0.129342553,  0.078231675,  0.081741171, -0.057731716, -0.025889917,
    0.009641678, -0.361916371,  0.138171131
  )
  
  
  
}

print("starting values:")
print(omega_choice_0)


# Full matrix to pass to the maximization routine:
data_max <- list(Nblocks = Nblocks, A = A, A_dum = A_dum, Achoice=Achoice, Z = Z, E=E, X = X, D = D, W = W,
                 A_big = A_big, Z_big = Z_big, pi_i = pi_i, lambda = lambda, homog = homog,
                 choice_col = choice_col, E_col = E_col, eta_out = eta_out,
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

# Maximization options (to be passed to optim)
optim_result <- optim(omega_choice_0,
                      fn = function(x) { objfn_no_partition(x, data_max, grad=FALSE)$Q },
                      gr = function(x) { objfn_no_partition(x, data_max, grad=TRUE)$G },
                      method = "L-BFGS-B", 
                      lower  = -Inf,
                      upper  =  Inf,
                      control = options_fmin)
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

if (homog == 0 && K > 1) {
  saveRDS(results_choice0_bfgs, file = "results_choice_bfgs.rds")
} else if (homog == 1) {
  saveRDS(results_choice0_bfgs, file = "results_choice_homog_bfgs.rds")
} else if (homog == 0 && K == 1) {
  saveRDS(results_choice0_bfgs, file = "results_choice_K1_bfgs.rds")
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
