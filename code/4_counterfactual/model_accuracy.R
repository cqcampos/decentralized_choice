################################################################################
# Script: run_counterfactual.R 
# Author: Chris Campos 
#
# Estimates counterfactual policy outcomes 
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
dir <- "/Volumes/lausd/decentralized_choice"
#dir <- "Z:/decentralized_choice"
#dir <- "/project/lausd/decentralized_choice"
setwd(dir)
FIGS <- file.path(dir, "output", "figures")

# Source helper files 
source(paste0(dir, "/code/helper/get_data.R"))
source(paste0(dir, "/code/helper/simulate_assignment_enrollment.R"))
source(paste0(dir, "/code/helper/calculate_eq.R"))

# Set up number of cores to use (adjust based on your machine)
#num_cores <- detectCores()

# Set parameters we will use in counterfactuals 
maxapps         <- 2
nsims           <- 200
preference_model <- 0   # 0=baseline, 1=lowcost, 2=altpref
accuracy       <- 1   # 1=Drop students used for MLE

# Set the model we are working with 
J <- 40
eta <- "eta_out"
homog <- "hetero"
K <- 3

################################################################################
############################ Read in data ######################################
################################################################################
# Load in seat capacity data 
seats <- read_dta(paste0(dir, "/data/prog_seats_2004_2008.dta"))
# Load structural sample
rawdata <- read_dta(paste0(dir, "/data/structural_data_2004_2008.dta"))

if (accuracy){
  # Load sample used for fitting MLE model 
  # 876 students submitted applications in multiple years (e.g., 2005 & 2006)
  # Among them, 302 students had one application used for MLE fitting and another not used
  # To assess model accuracy, we drop the applications used for MLE fitting
  # but retain the other applications, even if submitted by the same students
  rawdata_sample <- read_dta(paste0(dir, "/data/structural_data_2004_2008_sample.dta"))
  rawdata_sample$train <-TRUE
  rawdata <- left_join(rawdata, rawdata_sample[,c("studentpseudoid", "endyear", "train")], by = c("studentpseudoid", "endyear"))
  rawdata$train <- ifelse(is.na(rawdata$train), 0, rawdata$train)
  test_idx <- rawdata$train==FALSE
  rawdata <- rawdata[test_idx,]
  
  rawdata_sample <- NULL 
  gc()

}

# Load outcome model estimates
outcome_fes_math <- read.csv(paste0(dir, "/estimates/outcome_model_fes_",
                              eta, "_", homog, "_K", K, "_math.csv"))     # contains outcome fixed effects
outcome_ates_math <- read.csv(paste0(dir, "/estimates/outcome_model_ates_",
                              eta, "_", homog, "_K", K, "_math.csv"))   # contains outcome ates
outcome_betas_math <- read.csv(paste0(dir, "/estimates/outcome_model_betas_",
                              eta, "_", homog, "_K", K, "_math.csv"))   # contains outcome betas
outcome_fes_ela <- read.csv(paste0(dir, "/estimates/outcome_model_fes_",
                                    eta, "_", homog, "_K", K, "_ela.csv"))     # contains outcome fixed effects
names(outcome_fes_ela) <- c("endyear", "studentpseudoid", "yearfe_ela", "blockfe_ela")
outcome_ates_ela <- read.csv(paste0(dir, "/estimates/outcome_model_ates_",
                                     eta, "_", homog, "_K", K, "_ela.csv"))   # contains outcome ates
outcome_betas_ela <- read.csv(paste0(dir, "/estimates/outcome_model_betas_",
                                      eta, "_", homog, "_K", K, "_ela.csv"))   # contains outcome betas
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
pi_i <- as.matrix(read_dta(paste0(dir, "/data/student_p_i_2004_2008.dta")))
pi_i <- pi_i[test_idx, 3:ncol(pi_i)]
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
                                 homog, "_K", K, ".csv"))
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
  delta_j <- demand_params$omega[1:40]
  delta_X <- demand_params$omega[41:52]
  psi_D <- demand_params$omega[53:66]
  c_X <- demand_params$omega[67:79]
  
  sd_eta <- demand_params$omega[80]
  sd_theta_1 <- demand_params$omega[81]
  sd_theta_2 <- demand_params$omega[82]
  sd_theta_3 <- demand_params$omega[83]
  
  mu_theta_1 <- demand_params$omega[86]
  mu_theta_2 <- demand_params$omega[87]
  mu_theta_3 <- demand_params$omega[90]
  
  p_type_1 <- demand_params$omega[88]
  p_type_2 <- demand_params$omega[89]
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
seats <- read_dta(paste0(dir, "/data/prog_seats_2004_2008.dta"))

# Load in raw data
pi_i_init <- rawdata[, c("studentpseudoid", "endyear")]
pi_i_init$phbao <- (rawdata$white !=1)*1 # fed into calculations to merge on admission probs 


################################################################################
################## Simulate counterfactual outcomes  ###########################
################################################################################
results <- NULL
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
  pi_i_init
)


# 5.2 Repeat Monte Carlo draws
start_time <- Sys.time()
# Baseline 
params <- list(
  Amat            = Amat,
  pi_i            = pi_i,  # update later to eqm_pi[[j]],
  K               = K,
  homog           = homog,
  preference_model= 0
)

set.seed(123456)


# Helper function for constructing column idx of applied and enrolled school
find_col_index <- function(mat) {
  # Find column index of maximum value in each row
  result <- max.col(mat, ties.method = "first")
  # Handle rows where all values are 0
  # max.col will return 1 for all-zero rows, so we need to fix this
  row_sums <- rowSums(mat)
  result[row_sums == 0] <- -9999  
  return(result)
}

# Column index of program students applied to and predicted to apply to
A_expanded <- cbind((1-rowSums(A)), A)
choice_col <- find_col_index(A_expanded)
enr_col <- find_col_index(E)
applied <- (choice_col != 1)


# Create a placeholder for storing the simulation results

n_apps_pred <- matrix(nrow = nsims, ncol = (J+1)) # E(Number of applications to j)
n_offers_pred <- matrix(nrow = nsims, ncol = J) # E(Number of offers from j)
n_enr_pred <- matrix(nrow = nsims, ncol = J) # E(Number of students enrolling into j)

n_apps_phbao_pred <- matrix(nrow = nsims, ncol = (J+1)) # E(Number of applications to j) by phbao
n_offers_phbao_pred <- matrix(nrow = nsims, ncol = J) # E(Number of offers from j) by phbao
n_enr_phbao_pred <- matrix(nrow = nsims, ncol = J) # E(Number of students enrolling into j) by phbao

n_apps_white_pred <- matrix(nrow = nsims, ncol = (J+1)) # E(Number of applications to j) by white
n_offers_white_pred <- matrix(nrow = nsims, ncol = J) # E(Number of offers from j) by white
n_enr_white_pred <- matrix(nrow = nsims, ncol = J) # E(Number of students enrolling into j) by white

n_apps_asian_pred <- matrix(nrow = nsims, ncol = (J+1)) # E(Number of applications to j) by asian
n_offers_asian_pred <- matrix(nrow = nsims, ncol = J) # E(Number of offers from j) by asian
n_enr_asian_pred <- matrix(nrow = nsims, ncol = J) # E(Number of students enrolling into j) by asian

n_apps_low_score_school_pred <- matrix(nrow = nsims, ncol = (J+1)) # E(Number of applications to j) by low_score_school
n_offers_low_score_school_pred <- matrix(nrow = nsims, ncol = J) # E(Number of offers from j) by low_score_school
n_enr_low_score_school_pred <- matrix(nrow = nsims, ncol = J) # E(Number of students enrolling into j) by low_score_school

n_apps_high_score_school_pred <- matrix(nrow = nsims, ncol = (J+1)) # E(Number of applications to j) by high_score_school
n_offers_high_score_school_pred <- matrix(nrow = nsims, ncol = J) # E(Number of offers from j) by high_score_school
n_enr_high_score_school_pred <- matrix(nrow = nsims, ncol = J) # E(Number of students enrolling into j) by high_score_school

for (c in 1:nsims) {
  print(paste("Simulation", c, " of Baseline"))
  
  res <-  sim_a_e(
    endyear,
    low_score_school, 
    mean_school_scores,
    phbao, asian, white,
    yearfe,blockfe, outcome_fes, 
    alpha_j_math, alpha_0_X_math, alpha_m_X_math, alpha_0_math, alpha_m_math,
    alpha_j_ela, alpha_0_X_ela, alpha_m_X_ela, alpha_0_ela, alpha_m_ela,
    delta_j, delta_X, psi_D, c_X, 
    sd_eta, sd_theta_1, sd_theta_2, sd_theta_3, mu_theta_1, mu_theta_2, mu_theta_3,
    p_type_1, p_type_2,p_type_3,
    params,
    exog, eta 
  )
  
  n_apps_pred[c, ] <- res$n_apps
  n_offers_pred[c, ] <- res$n_offers
  n_enr_pred[c, ] <- res$n_enroll
  
  n_apps_phbao_pred[c,] <- res$n_apps_phbao 
  n_offers_phbao_pred[c,] <- res$n_offers_phbao
  n_enr_phbao_pred[c,] <- res$n_enroll_phbao
  
  n_apps_white_pred[c,] <- res$n_apps_white
  n_offers_white_pred[c,] <- res$n_offers_white
  n_enr_white_pred[c,] <- res$n_enroll_white
  
  n_apps_asian_pred[c,] <- res$n_apps_asian
  n_offers_asian_pred[c,] <- res$n_offers_asian
  n_enr_asian_pred[c,] <- res$n_enroll_asian
  
  n_apps_low_score_school_pred[c,] <- res$n_apps_low_score_school
  n_offers_low_score_school_pred[c,] <- res$n_offers_low_score_school
  n_enr_low_score_school_pred[c,] <- res$n_enroll_low_score_school
  
  n_apps_high_score_school_pred[c,] <- res$n_apps_high_score_school
  n_offers_high_score_school_pred[c,] <- res$n_offers_high_score_school
  n_enr_high_score_school_pred[c,] <- res$n_enroll_high_score_school
  
}

# Take the average over nsims

n_apps_avg <- colMeans(n_apps_pred)
n_offers_avg <- colMeans(n_offers_pred)
n_enr_avg    <- colMeans(n_enr_pred)


# Convert the numbers into shares
app_share_pred <- n_apps_avg/sum(n_apps_avg)
offer_rate_pred <- n_offers_avg/n_apps_avg[2:41]
enr_cond_offer_pred <- n_enr_avg/n_offers_avg

# Create sample analogues
app_share  <- colMeans(A_expanded)
offer_rate <- colSums(Z)/colSums(A)
enr_cond_offer <- colSums(E[,2:(J+1)])/colSums(Z)


app_shares <- data.frame(
  program = 0:J,
  predicted = app_share_pred,
  actual = app_share
)
app_shares <- app_shares[2:(J+1),]
write.csv(app_shares, file.path(dir, "estimates", paste0(eta, "_k", K, "_app_shares.csv")), row.names = FALSE)

offer_rates <- data.frame(
  program = 1:J,
  predicted = offer_rate_pred,
  actual = offer_rate
)
write.csv(offer_rates, file.path(dir, "estimates", paste0(eta, "_k", K, "_offer_rates.csv")), row.names = FALSE)

enr_rates <- data.frame(
  program = 1:J,
  predicted = enr_cond_offer_pred,
  actual = enr_cond_offer
)
write.csv(enr_rates, file.path(dir, "estimates", paste0(eta, "_k", K, "_enr_rates.csv")), row.names = FALSE)


# Create a csv file that combines all of the group-specific predictions and actuals 
# include apps, enrollment, offers
group_stats <- data.frame(
  program = 0:J,
  n_apps_pred = n_apps_avg,
  n_apps_actual = colSums(A_expanded),
  n_apps_phbao_pred = colMeans(n_apps_phbao_pred),
  n_apps_phbao_actual = colSums(A_expanded[which(phbao==1),]),
  n_apps_white_pred = colMeans(n_apps_white_pred),
  n_apps_white_actual = colSums(A_expanded[which(white==1),]),
  n_apps_asian_pred = colMeans(n_apps_asian_pred),
  n_apps_asian_actual = colSums(A_expanded[which(asian==1),]),
  n_apps_low_score_school_pred = colMeans(n_apps_low_score_school_pred),
  n_apps_low_score_school_actual = colSums(A_expanded[which(low_score_school==1),]),
  n_apps_high_score_school_pred = colMeans(n_apps_high_score_school_pred),
  n_apps_high_score_school_actual = colSums(A_expanded[which(low_score_school==0),])
)
group_stats2 <- data.frame(
  program = 1:J,
  n_offers_pred = n_offers_avg,
  n_offers_actual = colSums(Z),
  n_offers_phbao_pred = colMeans(n_offers_phbao_pred),
  n_offers_phbao_actual = colSums(Z[which(phbao==1),]),
  n_offers_white_pred = colMeans(n_offers_white_pred),
  n_offers_white_actual = colSums(Z[which(white==1),]),
  n_offers_asian_pred = colMeans(n_offers_asian_pred),
  n_offers_asian_actual = colSums(Z[which(asian==1),]),
  n_offers_low_score_school_pred = colMeans(n_offers_low_score_school_pred),
  n_offers_low_score_school_actual = colSums(Z[which(low_score_school==1),]),
  n_offers_high_score_school_pred = colMeans(n_offers_high_score_school_pred),
  n_offers_high_score_school_actual = colSums(Z[which(low_score_school==0),]),
  n_enr_pred = n_enr_avg,
  n_enr_actual = colSums(E[,2:(J+1)]),
  n_enr_phbao_pred = colMeans(n_enr_phbao_pred),
  n_enr_phbao_actual = colSums(E[which(phbao==1),2:(J+1)]),
  n_enr_white_pred = colMeans(n_enr_white_pred),
  n_enr_white_actual = colSums(E[which(white==1),2:(J+1)]),
  n_enr_asian_pred = colMeans(n_enr_asian_pred),
  n_enr_asian_actual = colSums(E[which(asian==1),2:(J+1)]),
  n_enr_low_score_school_pred = colMeans(n_enr_low_score_school_pred),
  n_enr_low_score_school_actual = colSums(E[which(low_score_school==1),2:(J+1)]),
  n_enr_high_score_school_pred = colMeans(n_enr_high_score_school_pred),
  n_enr_high_score_school_actual = colSums(E[which(low_score_school==0),2:(J+1)])
)
write.csv(group_stats, file.path(dir, "estimates", paste0(eta, "_k", K, "_group_stats.csv")), row.names = FALSE)
write.csv(group_stats2, file.path(dir, "estimates", paste0(eta, "_k", K, "_group_stats2.csv")), row.names = FALSE)
# Report averages from counterfactual and the sample analogue



bar_plot_vars <- function(v1, v2, lab1, lab2, xlab, ylab, out_path){

  
  # Combine into a matrix: each row is a vector, each column is an index
  mat <- rbind(v1, v2)
  
  pdf(out_path)
  
  # Create the bar plot
  barplot(mat,
          beside = TRUE,                # Put bars for the same index side-by-side
          col = c("skyblue", "orange"), # Colors for each vector
          names.arg = ((1:length(v1))-1), # X-axis labels
          legend.text = c(lab1, lab2), 
          xlab = xlab,
          ylab = ylab
          )
  
  dev.off()

}

# Portfolio level pred vs. actual
bar_plot_vars(app_share_pred, colMeans(A_expanded), "Predicted", "Actual", "Portfolio", "Share Choosing Portfolio", 
              file.path(FIGS, paste0(eta, "_k", K, "_share_bar.pdf")))

pdf(file.path(FIGS,paste0(eta, "_k", K, "_share_splot.pdf")))
plot(x=app_share_pred[2:41], y=colMeans(A),  xlab="Predicted Share", ylab="Actual Share")
dev.off()


pdf(file.path(FIGS,paste0(eta, "_k", K, "_share_splot_among_applicants.pdf")))
plot(x=n_apps_avg[2:41]/sum(n_apps_avg[2:41]),y=colMeans(A[applied,]),  xlab="Predicted Share", ylab="Actual Share")
dev.off()



# Enrollment | Predicted to get offer vs. Actual | Offer
pdf(file.path(FIGS, paste0(eta, "_k", K, "_share_splot_delta.pdf")))
plot(x=delta_j, y=colMeans(A), 
     xlab="j-specific mean utility", ylab="Actual Share")
dev.off()

pdf(file.path(FIGS, paste0(eta, "_k", K, "_pred_share_splot_delta.pdf")))
plot(x=delta_j, y=app_share_pred[2:41], 
     xlab="j-specific mean utility", ylab="Predicted Share")
dev.off()


# Offer rate, conditional on applying
pdf(file.path(FIGS, paste0(eta, "_k", K, "_offer_rate.pdf")))
plot(y=offer_rate, x=offer_rate_pred, 
     ylab="Actual Offer Rate | Applying", 
     xlab="Predicted Offer Rate | Applying")
dev.off()

# Enrollment rate, conditional on offer
pdf(file.path(FIGS, paste0(eta, "_k", K, "_e_rate.pdf")))
plot(y=enr_cond_offer, x=enr_cond_offer_pred,
     ylab="Actual Enrollment Rate | Offer",
     xlab="Predicted Enrollment Rate | Offer")
dev.off()
