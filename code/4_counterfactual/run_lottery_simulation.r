################################################################################
# Script: run_lottery_simulation.r 
# Author: Chris Campos
# Last Edited: Chris Campos, October 2025 
# Estimates counterfactual lotteries  
# 
# Dependencies:

################################################################################
rm(list = ls())
gc()
options(error = recover) 



################################################################################
########################### Preliminaries ######################################
################################################################################
# Set personal library path
#install.packages(c("haven", "doParallel", "foreach", "numDeriv", "parallel"))

# Establish paths 
dir <- "/Volumes/lausd/decentralized_choice/code/chris/decentralized_choice"
# dir <- "Z:/decentralized_choice"
#dir <- "/project/lausd/decentralized_choice"
win_os <- (dir =="Z:/decentralized_choice")


library(parallel)    # For parallel processing
library(numDeriv)    # For numerical derivatives
library(haven)       # For reading Stata files
library(doParallel)  # For foreach parallel backend
library(foreach)     # For parallel loops
library(dplyr)
library(tidyr)
library(matrixStats) # row wise, column wise stats
library(doParallel) # For checking number of cores available
library(reticulate) # Python


# Source helper files 
source(paste0(dir, "/code/helper/get_data.R"))



# Set parameters we will use in counterfactuals 
maxapps         <- 1
nsims           <- 100

# Set the model we are working with 
J <- 53
last_yr <- 2013
app_2013 <- 1
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
rawdata$low_score_school <- ifelse(rawdata$mean_school_scores < -0.3, 1, 0)
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
alpha_0_X_math <- outcome_betas_math$est_math[1:11]
alpha_m_X_math <- outcome_betas_math$est_math[12:23]
alpha_0_math <- outcome_betas_math$est_math[24]
alpha_m_math <- outcome_betas_math$est_math[25]
# Outcome model: ela
alpha_j_ela <- outcome_ates_ela$estimate
alpha_0_X_ela <- outcome_betas_ela$est_ela[1:11]
alpha_m_X_ela <- outcome_betas_ela$est_ela[12:23]
alpha_0_ela <- outcome_betas_ela$est_ela[24]
alpha_m_ela <- outcome_betas_ela$est_ela[25]

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
  endyear = endyear,
  pi_i = pi_i,
  seats_df = seats,
  phbao = phbao,
  eta             = eta,
  param_init       = param,
  last_yr          = last_yr
)

# params (you supply; names match the script)
params <- list(
  # utilities
  delta_j  = delta_j,              
  delta_X  = delta_X,             
  psi_D    = psi_D,                 

  # app costs
  c_X      = c_X,                   
  sd_eta   = sd_eta,                 

  # types
  K        = 3,
  mu_theta = c(mu_theta_1, mu_theta_2, mu_theta_3),
  sd_theta = c(sd_theta_1, sd_theta_2, sd_theta_3),    
  p_type   = c(p_type_1, p_type_2, p_type_3),

  # outcomes
  alpha_j_math = alpha_j_math, 
  alpha_0_X_math = alpha_0_X_math, 
  alpha_m_X_math = alpha_m_X_math, 
  alpha_0_math = alpha_0_math, 
  alpha_m_math = alpha_m_math,           
  alpha_j_ela = alpha_j_ela, 
  alpha_0_X_ela = alpha_0_X_ela, 
  alpha_m_X_ela = alpha_m_X_ela, 
  alpha_0_ela = alpha_0_ela, 
  alpha_m_ela = alpha_m_ela,
  sd_eps0    = 1.0,
  sd_eps1    = 1.0
)

source(paste0(dir, "/code/4_counterfactual/simulate_lotteries.R"))

# Repeat Monte Carlo draws
panel <- simulate_lottery_panel(
  nsims   = 1,
  exog    = exog,
  params  = params,
  out_file = paste0(dir, "/estimates/sim_lottery_panel.csv"),
  restrict_to_oversubscribed = FALSE,
  scale_app = 0.05,
  seed = 12345
)
