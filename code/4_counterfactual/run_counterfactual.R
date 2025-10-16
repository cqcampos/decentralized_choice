################################################################################
# Script: run_counterfactual.R 
# Author: Chris Campos, Ryan Lee
# Last Edited: Ryan, October 16, 2025 
# Estimates counterfactual policy outcomes 
# 
# Dependencies:
################################################################################
rm(list = ls())
gc()
options(error = recover) 

################################################################################
########################### Preliminaries ######################################
################################################################################
# --- Establish paths  ---
slurm_user <- Sys.getenv("SLURM_JOB_USER") 
if (slurm_user == "") {
  # When running code in local device ---
  # Detect OS and set base directory
  os_name <- Sys.info()[["sysname"]]
  if (os_name == "Windows") {
    dir <- "Z:/decentralized_choice"
  } else {
    dir <- "/Volumes/lausd/decentralized_choice"
  }
} else {
  # When running code via Slurm ---
  # Update path based on Slurm account
  dir <- "/project/lausd/decentralized_choice"
  
  # Update cloned repo path to run code from
  if (slurm_user == "faculty") {
    code_dir <- "/project/lausd/decentralized_choice/code/chris/decentralized_choice"
  } else if (slurm_user == "ryanlee22") {
    code_dir <- "/project/lausd/decentralized_choice/code/ryan/decentralized_choice"
  } else {
    stop(sprintf(
      "Unrecognized Slurm user: '%s'.\nPlease update `code_dir` to point to your cloned repository on the server.",
      slurm_user
    ))
  }
}
# For determining method of parallel computation
win_os <- (dir =="Z:/decentralized_choice")

# Use torch for finding equilibrium pi
use_torch_for_eq <- TRUE

# Use the mean(pi_br), instead of optimizing each round
use_mean_pi      <-  FALSE 


# --- Load library and helper functions ---
# Point to the Python binary from the module system
# Failing to specify python path can cause issue when reading NPZ file:
# reticulate won't convert the numpy array to matrix. 
use_python("/apps/python/3.10/3.10.9/bin/python3", required = TRUE) # For Mercury Cluster Python 3.10
reticulate::py_install("numpy")
reticulate::py_install("numba")

# When running code in your local device:
# reticulate::install_python("3.12")
# use_python("/usr/bin/python3", required = TRUE) # Yes, for mac users 
# use_python("C:/Users/Administrator/AppData/Local/Programs/Python/Python312/python.exe", required = TRUE) # For windows
# use_python("C:/Users/{USERNAME}/AppData/Local/r-reticulate/r-reticulate/pyenv/pyenv-win/versions/3.12.3/python.exe", required = TRUE) # For windows


# Load R libraries
if (use_torch_for_eq){
  library("torch")
  print("Is CUDA available?")
  print(cuda_is_available())
}
library(parallel)    # For parallel processing
library(numDeriv)    # For numerical derivatives
library(haven)       # For reading Stata files
library(doParallel)  # For foreach parallel backend
library(foreach)     # For parallel loops
library(dplyr)
library(tidyr)
library(matrixStats) # row wise, column wise stats
library(doParallel)  # For checking number of cores available
library(reticulate)  # Python


# Source helper files 
source(paste0(code_dir, "/code/helper/get_data.R"))
source(paste0(code_dir,"/code/helper/simulate_counterfactual.r")) 
source(paste0(code_dir, "/code/helper/pi_norm_torch.R"))
source(paste0(code_dir, "/code/helper/calculate_eq.R"))
source(paste0(code_dir, "/code/helper/matching.R"))


# --- Log file setup ---
log_file <- paste0(code_dir,  "/code/4_counterfactual/logs/run_counterfactual_", format(Sys.time(), "%Y_%m_%d_%H_%M_%S"), ".txt")
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
log_message("Starting counterfactual estimation...")

# --- Set parameters we will use in counterfactuals ---
maxapps         <- 1
nsims           <- 100

# Sample selection
accuracy       <- 0   # 1=Drop students used for MLE
app_2013       <- 1   # 1=use 2004-13 sample, 0=2004-08 sample

# Set the model we are working with 
J <- ifelse(app_2013, 53, 40)
last_yr <- ifelse(app_2013, 2013, 2008)
eta <- "eta_out"
homog <- "hetero"
K <- 3

print(paste("Model:", eta))
print(paste("K:", K))

################################################################################
############################ Read in data ######################################
################################################################################
# Load in seat capacity data 
seats <- read_dta(paste0(dir, "/data/prog_seats_2004_", last_yr, ".dta"))
# Load structural sample
rawdata <- read_dta(paste0(dir, "/data/structural_data_2004_", last_yr, ".dta"))
rawdata$districtdum2 <- ifelse(rawdata$localdistrictcode=="E", 1, 0)
rawdata$districtdum3 <- ifelse(rawdata$localdistrictcode=="NE",1, 0)
rawdata$districtdum4 <- ifelse(rawdata$localdistrictcode=="NW",1, 0)
rawdata$districtdum5 <- ifelse(rawdata$localdistrictcode=="S",1, 0)
rawdata$districtdum6 <- ifelse(rawdata$localdistrictcode=="W",1, 0)

if (accuracy){
  # Load sample used for fitting MLE model 
  # 876 students submitted applications in multiple years (e.g., 2005 & 2006)
  # Among them, 302 students had one application used for MLE fitting and another not used
  # To assess model accuracy, we drop the applications used for MLE fitting
  # but retain the other applications, even if submitted by the same students
  rawdata_sample <- read_dta(paste0(dir, "/data/structural_data_2004_", last_yr, "_sample.dta"))
  rawdata_sample$train <-TRUE
  rawdata <- left_join(rawdata, rawdata_sample[,c("studentpseudoid", "endyear", "train")], by = c("studentpseudoid", "endyear"))
  rawdata$train <- ifelse(is.na(rawdata$train), 0, rawdata$train)
  rawdata <- rawdata[rawdata$train==FALSE,]
  
  rawdata_sample <- NULL 
  gc()

}

# Load outcome model estimates
app_suffix <- ifelse(app_2013, "_2013", "")
outcome_fes_math <- read.csv(paste0(dir, "/estimates/outcome_model_fes_",
                              eta, "_", homog, "_K", K, app_suffix, "_math.csv"))     # contains outcome fixed effects
outcome_ates_math <- read.csv(paste0(dir, "/estimates/outcome_model_ates_",
                              eta, "_", homog, "_K", K, app_suffix, "_math.csv"))   # contains outcome ates
outcome_betas_math <- read.csv(paste0(dir, "/estimates/outcome_model_betas_",
                              eta, "_", homog, "_K", K, app_suffix, "_math.csv"))   # contains outcome betas
outcome_fes_ela <- read.csv(paste0(dir, "/estimates/outcome_model_fes_",
                                    eta, "_", homog, "_K", K, app_suffix, "_ela.csv"))     # contains outcome fixed effects
names(outcome_fes_ela) <- c("endyear", "studentpseudoid", "yearfe_ela", "blockfe_ela")
outcome_ates_ela <- read.csv(paste0(dir, "/estimates/outcome_model_ates_",
                                     eta, "_", homog, "_K", K, app_suffix, "_ela.csv"))   # contains outcome ates
outcome_betas_ela <- read.csv(paste0(dir, "/estimates/outcome_model_betas_",
                                      eta, "_", homog, "_K", K, app_suffix, "_ela.csv"))   # contains outcome betas
# Merge outcome model fixed effects and create a few other useful variables 
outcome_fes_ela$studentpseudoid <- as.character(outcome_fes_ela$studentpseudoid)
outcome_fes_math$studentpseudoid <- as.character(outcome_fes_math$studentpseudoid)
rawdata <- rawdata %>%
  inner_join(outcome_fes_math, by = c("studentpseudoid" = "studentpseudoid", "endyear" = "endyear")) %>%
  inner_join(outcome_fes_ela, by = c("studentpseudoid" = "studentpseudoid", "endyear" = "endyear"))
rawdata <- rawdata %>% group_by(preferredlocationcode, endyear) %>%
  mutate(mean_school_scores = mean(lag_math, na.rm=TRUE))
rawdata$low_score_school <- ifelse(rawdata$mean_school_scores < -0, 1, 0)
white <- rawdata$white
asian <- rawdata$asian
phbao <- 1- rawdata$white
yearfe <- rawdata$yearfe
blockfe <- rawdata$blockfe
yearfe_ela <- rawdata$yearfe_ela
blockfe_ela <- rawdata$blockfe_ela
asian <- rawdata$asian
endyear <- rawdata$endyear
mean_school_scores <- rawdata$mean_school_scores
low_score_school <- rawdata$low_score_school

# Use helper function to create matrices relevant for counterfactuals 
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


# Load admission probabilities and the admission matrix:
pi_i <- as.matrix(read_dta(paste0(dir, "/data/student_p_i_2004_", last_yr, ".dta")))
pi_i <- pi_i[, 3:ncol(pi_i)]
storage.mode(pi_i) <- "numeric"

# Additional prep on some matrices 
X <- scale(X, center = TRUE, scale = FALSE)
Q_X <- ncol(X)

# Additional variables used for the outcome model
X_Y <- as.matrix(cbind( rawdata$missing_ela, rawdata$missing_math,
                       rawdata$districtdum2, 
                       rawdata$districtdum3,
                       rawdata$districtdum4,
                       rawdata$districtdum5,
                       rawdata$districtdum6))
X_Y <- scale(X_Y, center = TRUE, scale = FALSE)

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

# Compute dims
N   <- nrow(Y)
Q_X <- ncol(X)
# D is N × Q_D × J, so Q_D = dim(D)[2]
Q_D <- dim(D)[2]
Q_W <- ncol(W)


################################################################################
##################### Estimated Model Parameters ###############################
################################################################################
# Demand model 
demand_params <- read.csv(paste0(dir, "/estimates/in_mag_ind_", 
                                 eta, "_estimates_choice_",
                                 homog, "_K", K, app_suffix, ".csv"))
if(K==1){
  delta_j <- demand_params$omega[1:40]
  delta_X <- demand_params$omega[41:52]
  psi_D <- demand_params$omega[53:66]
  c_X <- demand_params$omega[69:81]
  sd_eta <- demand_params$omega[67]
  sd_theta_1 <- demand_params$omega[68]
  sd_theta_2 <- 0 
  sd_theta_3 <- 0 
  mu_theta_1 <- 0
  mu_theta_2 <- 0
  mu_theta_3 <- 0
  p_type_1 <- 1
  p_type_2 <- 0
  p_type_3 <- 0
}
if(K==2){
  delta_j <- demand_params$omega[1:40]
  delta_X <- demand_params$omega[41:52]
  psi_D <- demand_params$omega[53:66]
  c_X <- demand_params$omega[67:79]
  
  sd_eta <- demand_params$omega[80]
  sd_theta_1 <- demand_params$omega[81]
  sd_theta_2 <- demand_params$omega[82]
  sd_theta_3 <- 0
  
  mu_theta_1 <- demand_params$omega[84]
  mu_theta_2 <- demand_params$omega[86]
  mu_theta_3 <- 0
  
  p_type_1 <- demand_params$omega[85]
  p_type_2 <- 1-p_type_1
  p_type_3 <- 0 
}
if(K==3){
  # 13 more magnet schools when we use 2004-2013 sample
  a <- ifelse(app_2013, 13, 0)
  delta_j <- demand_params$omega[1:(40+a)]
  delta_X <- demand_params$omega[(41+a):(52+a)]
  psi_D <- demand_params$omega[(53+a):(66+a)]
  c_X <- demand_params$omega[(67+a):(79+a)]
  
  sd_eta <- demand_params$omega[80+a]
  sd_theta_1 <- demand_params$omega[81+a]
  sd_theta_2 <- demand_params$omega[82+a]
  sd_theta_3 <- demand_params$omega[83+a]
  
  mu_theta_1 <- demand_params$omega[86+a]
  mu_theta_2 <- demand_params$omega[87+a]
  mu_theta_3 <- demand_params$omega[90+a]
  
  p_type_1 <- demand_params$omega[88+a]
  p_type_2 <- demand_params$omega[89+a]
  p_type_3 <- 1 - p_type_1 - p_type_2
}
# Outcome model: math
alpha_j_math <- outcome_ates_math$estimate 
alpha_0_X_math <- outcome_betas_math$est_math[1:13]
alpha_m_X_math <- outcome_betas_math$est_math[14:30]
alpha_0_math <- outcome_betas_math$est_math[32]
alpha_m_math <- outcome_betas_math$est_math[31]
# Outcome model: ela
alpha_j_ela <- outcome_ates_ela$estimate
alpha_0_X_ela <- outcome_betas_ela$est_ela[1:13]
alpha_m_X_ela <- outcome_betas_ela$est_ela[14:30]
alpha_0_ela <- outcome_betas_ela$est_ela[32]
alpha_m_ela <- outcome_betas_ela$est_ela[31]

outcome_params <- list(
  alpha_j_math = alpha_j_math,
  alpha_0_X_math = alpha_0_X_math,
  alpha_m_X_math = alpha_m_X_math,
  alpha_0_math = alpha_0_math,
  alpha_m_math = alpha_m_math,
  alpha_j_ela = alpha_j_ela,
  alpha_0_X_ela = alpha_0_X_ela,
  alpha_m_X_ela = alpha_m_X_ela,
  alpha_0_ela = alpha_0_ela,
  alpha_m_ela = alpha_m_ela
)

################################################################################
################## Prep eq. admission prob outcomes  ###########################
################################################################################
seats_eq <- read_dta(paste0(dir, "/data/prog_seats_2004_", last_yr, ".dta"))
# Reshape so that we have one row per school per year per phbao
seats_eq <- seats_eq %>%
  pivot_longer(
    cols      = starts_with("seats_"),     # all seat cols
    names_to  = "schoolcode",               # new key column
    names_prefix = "seats_",               # drop the prefix
    values_to = "seats"                    # new value column
  ) %>% group_by(schoolcode, endyear) %>%
  mutate(total_seats = sum(seats)) %>%
  arrange(schoolcode, endyear)
# Calculate share of seats distributed to phbao/others for each school x endyear
seats_eq$share_seats <- seats_eq$seats / seats_eq$total_seats

# Create a student x endyear level datasets telling whether a student is phbao
stu_type <- rawdata[, c("studentpseudoid", "endyear")]
stu_type$phbao <- (rawdata$white !=1)*1 # fed into calculations to merge on admission probs 
stu_type <- as.data.frame(stu_type)

################################################################################
################## Simulate counterfactual outcomes  ###########################
################################################################################
# Save starting values for admission prob eq. 
param <- read_dta(paste0(dir, "/data/aggregate_prog_admission_prob_2004_", last_yr, ".dta"))

# Estimate 1 per school x endyear x phbao
app_years <- sort(unique(stu_type$endyear))
n_types <- length(app_years) * 2          
param <- rep(as.numeric(param$p), each = n_types)

# Matrices of possible A and Z values
# Can only receive one offer 
Zmat <- diag(J)
Zmat <- rbind(Zmat, rep(0, J))
Amat <- Zmat
Q_A <- nrow(Amat)
Q_Z <- nrow(Zmat)
AZ <- list(Amat = Amat, Zmat = Zmat)


exog <- list(
  X = X,
  D = D,
  W = W,
  seats=seats, 
  seats_eq = seats_eq, 
  stu_type =stu_type,
  eta             = eta,
  param_init       = param,
  maxapps          = maxapps,
  last_yr          = last_yr
)

# Repeat Monte Carlo draws
# p_mod == 0 <- baseline
# p_mod == 1 <- No application cost
# p_mod == 2 <- Info provision
# p_mod == 3 <- Targeted info provision
# p_mod == 4 <- Unified enrollments
# p_mod == 5 <- remove travel costs 
# p_mod == 6 <- DA + no travel costs  
# p_mod == 7 <- optimal sorting on match quality 
# p_mod == 8 <- Unified Enrollment (DA) + Info provision 
# p_mod == 9 <- Unified Enrollment (DA) + no travel costs + Info Provision
for (p_mod in c(0)){
  run_sim_cf(dir, win_os, nsims, p_mod, endyear, low_score_school, mean_school_scores, 
             phbao, asian, white, 
             yearfe, blockfe, outcome_fes, 
             alpha_j_math, alpha_0_X_math, alpha_m_X_math, alpha_0_math, alpha_m_math,
             alpha_j_ela, alpha_0_X_ela, alpha_m_X_ela, alpha_0_ela, alpha_m_ela,
             delta_j, delta_X, psi_D, c_X,  sd_eta,
             sd_theta_1, sd_theta_2, sd_theta_3,
             mu_theta_1, mu_theta_2, mu_theta_3,
             p_type_1, p_type_2, p_type_3,
             params, exog, use_torch_for_eq, use_mean_pi =use_mean_pi)

}
