################################################################################
# Script: admission_prob_eq.R 
# Author: Chris Campos 
#
# Estimates equilibrium admission probs 
# 
# Dependencies:

################################################################################
rm(list = ls())
gc()

################################################################################
########################### Preliminaries ######################################
################################################################################
# Set personal library path
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
#install.packages(c("haven", "doParallel", "foreach", "numDeriv", "parallel"))

# Load necessary libraries
library(parallel)    # For parallel processing
library(numDeriv)    # For numerical derivatives
library(haven)       # For reading Stata files
library(doParallel)  # For foreach parallel backend
library(foreach)     # For parallel loops
library(dplyr)
library(tidyr)
library(matrixStats)

# Establish paths 
setwd("/project/lausd/decentralized_choice")
dir <- "/project/lausd/decentralized_choice"
# Set the model we are working with 
J <- 40
eta <- "eta_in"
homog <- "hetero"
K <- 2
R <- 200

# Source helper files 
source(paste0(dir, "/code/get_data.R"))
source(paste0(dir, "/code/admission_prob_norm.R"))

################################################################################
########################### Prepare Data #######################################
################################################################################
# Load in raw data
rawdata <- read_dta(paste0(dir, "/data/structural_data_2004_2008.dta"))
N <- nrow(rawdata)
pi_i <- rawdata[, c("studentpseudoid", "endyear")]
pi_i$phbao <- (rawdata$white !=1)*1 # fed into calculations to merge on admission probs 

# Use helper function to create matrices relevant for calculations 
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




# Load in seat capacity data and prepare for admission prob eq calculations 
seats <- read_dta(paste0(dir, "/data/prog_seats_2004_2008.dta"))
# Reshape so that we have one row per school per year per phbao
seats_long <- seats %>%
  pivot_longer(
    cols      = starts_with("seats_"),     # all 40 seat cols
    names_to  = "schoolcode",               # new key column
    names_prefix = "seats_",               # drop the prefix
    values_to = "seats"                    # new value column
  ) %>% group_by(schoolcode) %>%
  mutate(total_seats = sum(seats)) 
# Distribute seats within a school across years and phbao to determine seats 
seats_long$share_seats <- seats_long$seats / seats_long$total_seats
seats <- seats_long %>%
  # now group at the level you actually want …
  group_by(endyear, schoolcode, phbao) %>%  
  # … and sum up your shares (and, if you like, the raw seats too)
  summarize(
    total_seats    = sum(seats,       na.rm = TRUE),
    share_seats    = sum(share_seats, na.rm = TRUE)
  ) %>% arrange(schoolcode)
rm(seats_long)

################################################################################
##################### Estimated Model Parameters ###############################
################################################################################
# Demand model 
demand_params <- read.csv(paste0(dir, "/estimates/", 
                                 eta, "_estimates_choice_",
                                 homog, "_K", K, ".csv"))
delta_j <- demand_params$omega[1:40]
delta_X <- demand_params$omega[41:52]
psi_D <- demand_params$omega[53:66]
c_X <- demand_params$omega[67:79]

sd_eta <- demand_params$omega[80]
sd_theta_1 <- demand_params$omega[81]
sd_theta_2 <- demand_params$omega[82]
mu_theta_1 <- demand_params$omega[84]
mu_theta_2 <- demand_params$omega[86]
p_type_1 <- demand_params$omega[85]
p_type_2 <- 1-p_type_1

preference_model <-1
param <- read_dta(paste0(dir, "/data/aggregate_prog_admission_prob_2004_2008.dta"))
param <- param$p
param <- as.numeric(param)*3
#param <- rep(10, 40)
options_fmin <- list(maxit = 100000, reltol = 1e-8)

pi_br_optim <- optim(param,
                        fn = function(x) { admission_prob_norm(x, seats=seats,
                                                               pi_i_init = pi_i,
                                                               delta_j = delta_j, delta_X = delta_X,
                                                               psi_D = psi_D, c_X = c_X,
                                                               sd_eta = sd_eta,
                                                               sd_theta_1 = sd_theta_1,
                                                               sd_theta_2 = sd_theta_2,
                                                               mu_theta_1 = mu_theta_1,
                                                               mu_theta_2 = mu_theta_2,
                                                               p_type_1 = p_type_1,
                                                               p_type_2 = p_type_2,
                                                               preference_model = preference_model, R=R)
                                                        },
                        method = "Nelder-Mead", 
                        lower  = -Inf,
                        upper  =  Inf,
                        control = options_fmin)

print(pi_br_optim$par)


