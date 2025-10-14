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
setwd(dir)
getwd()
# Source helper files 
source(paste0(dir, "/code/get_data.R"))
source(paste0(dir,"/code/simulate_counterfactual.r")) 
source(paste0(dir, "/code/pi_norm.R"))
source(paste0(dir, "/code/calculate_eq.R"))
# Set up number of cores to use (adjust based on your machine)
#num_cores <- detectCores()

################################# Log file setup ###############################
log_file <- paste0(getwd(),  "/logs/counterfactuals_", format(Sys.time(), "%Y_%m_%d_%H_%M_%S"), ".txt")
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

# Set parameters we will use in counterfactuals 
maxapps         <- 2
nsims           <- 200
preference_model <- 0   # 0=baseline, 1=lowcost, 2=altpref
accuracy       <- 0   # 1=Drop students used for MLE

# Set the model we are working with 
J <- 40
eta <- "eta_in"
homog <- "hetero"
K <- 2

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
  rawdata <- rawdata[rawdata$train==FALSE,]
  
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
seats_eq <- read_dta(paste0(dir, "/data/prog_seats_2004_2008.dta"))
# Reshape so that we have one row per school per year per phbao
seats_long <- seats_eq %>%
  pivot_longer(
    cols      = starts_with("seats_"),     # all 40 seat cols
    names_to  = "schoolcode",               # new key column
    names_prefix = "seats_",               # drop the prefix
    values_to = "seats"                    # new value column
  ) %>% group_by(schoolcode) %>%
  mutate(total_seats = sum(seats)) 
# Distribute seats within a school across years and phbao to determine seats 
seats_long$share_seats <- seats_long$seats / seats_long$total_seats
seats_eq <- seats_long %>%
  # now group at the level you actually want …
  group_by(endyear, schoolcode, phbao) %>%  
  # … and sum up your shares (and, if you like, the raw seats too)
  summarize(
    total_seats    = sum(seats,       na.rm = TRUE),
    share_seats    = sum(share_seats, na.rm = TRUE)
  ) %>% arrange(schoolcode)
rm(seats_long)
# Load in raw data
pi_i_init <- rawdata[, c("studentpseudoid", "endyear")]
pi_i_init$phbao <- (rawdata$white !=1)*1 # fed into calculations to merge on admission probs 


################################################################################
################## Simulate counterfactual outcomes  ###########################
################################################################################

# Save starting values for admission prob eq. 
param <- read_dta(paste0(dir, "/data/aggregate_prog_admission_prob_2004_2008.dta"))
param <- param$p
param <- as.numeric(param)*3
write.csv(param, paste0(dir, "/estimates/init_param_preference_model_0.csv"), row.names = FALSE)
write.csv(param, paste0(dir, "/estimates/init_param_preference_model_1.csv"), row.names = FALSE)
write.csv(param, paste0(dir, "/estimates/init_param_preference_model_2.csv"), row.names = FALSE)
write.csv(param, paste0(dir, "/estimates/init_param_preference_model_3.csv"), row.names = FALSE)

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
  E = E,
  seats = seats, 
  seats_eq = seats_eq, 
  pi_i_init = pi_i_init
)

# FIXED: Create outcome_fes if it doesn't exist
outcome_fes <- list(
  yearfe = yearfe,
  blockfe = blockfe,
  yearfe_ela = yearfe_ela,
  blockfe_ela = blockfe_ela
)

library(doParallel)
library(foreach)
library(doRNG)

log_message("Starting Parallelization Prep...")

scenarios <- data.frame(
  preference_model = 0:3,
  label = c("Baseline", "NoAppCosts", "InfoLight", "InfoTargeted")
)
task_df <- merge(scenarios, data.frame(sim = seq_len(nsims)), all = TRUE)

fixed_params <- list(Amat = Amat, pi_i = pi_i, K = K, homog = homog)

# Set up cluster
ncores <- max(1L, parallel::detectCores() - 1L)
log_message(paste("Using", ncores, "cores for parallel processing"))
cl <- makeCluster(ncores, outfile = "") # outfile="" shows worker output
registerDoParallel(cl)

# Set RNG for reproducibility
RNGkind("L'Ecuyer-CMRG")
set.seed(123456)

# Create separate log files for each worker to avoid file locking
worker_logs_dir <- paste0(dir, "/logs/workers/")
if (!dir.exists(worker_logs_dir)) {
  dir.create(worker_logs_dir, recursive = TRUE)
}

# Export all required variables to workers
clusterExport(cl, c(
  # Data
  "endyear", "low_score_school", "mean_school_scores", 
  "phbao", "asian", "white", "yearfe", "blockfe", 
  "yearfe_ela", "blockfe_ela", "outcome_fes",
  # Parameters
  "alpha_j_math", "alpha_0_X_math", "alpha_m_X_math", "alpha_0_math", "alpha_m_math",
  "alpha_j_ela", "alpha_0_X_ela", "alpha_m_X_ela", "alpha_0_ela", "alpha_m_ela",
  "delta_j", "delta_X", "psi_D", "c_X",
  "sd_eta", "sd_theta_1", "sd_theta_2", "mu_theta_1", "mu_theta_2", 
  "p_type_1", "p_type_2",
  # Other
  "fixed_params", "exog", "dir", "worker_logs_dir"
))

# Function to monitor progress
monitor_progress <- function(task_df, worker_logs_dir, total_time = 60) {
  start_monitor <- Sys.time()
  completed <- rep(FALSE, nrow(task_df))
  
  while(sum(completed) < nrow(task_df)) {
    # Check for completed tasks by looking for result files
    for(i in 1:nrow(task_df)) {
      if(!completed[i]) {
        result_file <- paste0(worker_logs_dir, "task_", i, "_complete.txt")
        if(file.exists(result_file)) {
          completed[i] <- TRUE
          log_message(paste("Task", i, "completed:",
                          "Model =", task_df$label[i],
                          "Sim =", task_df$sim[i],
                          "| Progress:", sum(completed), "/", nrow(task_df)))
        }
      }
    }
    
    # Check timeout
    if(difftime(Sys.time(), start_monitor, units = "hours") > total_time) {
      log_message("WARNING: Monitoring timeout reached")
      break
    }
    
    Sys.sleep(10) # Check every 10 seconds
  }
}

# Start monitoring in background (optional - you can run this manually)
# monitor_job <- parallel::mcparallel(monitor_progress(task_df, worker_logs_dir))

log_message(paste("Starting", nrow(task_df), "tasks..."))
log_message("Baseline tasks (preference_model=0) should complete quickly")
log_message("Other tasks may take 15-60 minutes each")

# Run parallel tasks with better error handling
res_list <- foreach(
  k = seq_len(nrow(task_df)),
  .combine = 'list',  # Use list first, then combine
  .inorder = FALSE,
  .errorhandling = "pass",  # Return errors instead of removing
  .options.RNG = 123456,
  .export = c("simulate_counterfactual", "calculate_eq", "pi_norm"),
  .packages = c("dplyr", "matrixStats", "haven", "tidyr")
) %dorng% {
  
  # Create worker-specific log file
  worker_log_file <- paste0(worker_logs_dir, "worker_task_", k, ".log")
  
  # Simple logging function for this worker
  worker_log <- function(msg) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    formatted_msg <- paste0("[Task ", k, " - ", timestamp, "] ", msg)
    cat(formatted_msg, "\n", file = worker_log_file, append = TRUE)
  }
  
  tryCatch({
    start_time <- Sys.time()
    worker_log(paste0("Starting task: Model = ", task_df$label[k], 
                     ", Sim = ", task_df$sim[k]))
    
    pm <- task_df$preference_model[k]
    sim <- task_df$sim[k]
    params <- c(fixed_params, list(preference_model = pm))
    
    # Log before calling the main function
    worker_log("Calling simulate_counterfactual...")
    
    stats <- simulate_counterfactual(
      endyear,
      low_score_school, 
      mean_school_scores,
      phbao, asian, white,
      yearfe, blockfe, outcome_fes, 
      alpha_j_math, alpha_0_X_math, alpha_m_X_math, alpha_0_math, alpha_m_math,
      alpha_j_ela, alpha_0_X_ela, alpha_m_X_ela, alpha_0_ela, alpha_m_ela,
      delta_j, delta_X, psi_D, c_X, 
      sd_eta, sd_theta_1, sd_theta_2, mu_theta_1, mu_theta_2, p_type_1, p_type_2,
      params,
      exog
    )
    
    end_time <- Sys.time()
    duration <- difftime(end_time, start_time, units = "mins")
    worker_log(paste0("Completed in ", round(duration, 2), " minutes"))
    
    # Write completion marker
    writeLines("complete", paste0(worker_logs_dir, "task_", k, "_complete.txt"))
    
    # Return results
    list(
      success = TRUE,
      sim = sim,
      preference_model = pm,
      stats = stats,
      duration = as.numeric(duration)
    )
    
  }, error = function(e) {
    worker_log(paste0("ERROR: ", e$message))
    list(
      success = FALSE,
      sim = task_df$sim[k],
      preference_model = task_df$preference_model[k],
      error = e$message
    )
  })
}

stopCluster(cl)

# Process results
successful_results <- Filter(function(x) x$success, res_list)
failed_results <- Filter(function(x) !x$success, res_list)

log_message(paste("Completed:", length(successful_results), "successful,", 
                 length(failed_results), "failed"))

if(length(failed_results) > 0) {
  log_message("Failed tasks:")
  for(f in failed_results) {
    log_message(paste("  Sim", f$sim, "Model", f$preference_model, ":", f$error))
  }
}

# Convert successful results to matrix
if(length(successful_results) > 0) {
  res_mat <- do.call(rbind, lapply(successful_results, function(x) {
    c(sim = x$sim, preference_model = x$preference_model, x$stats)
  }))
  
  # Save results
  outfile <- paste0(dir, "/estimates/counterfactual_results2_", 
                   format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
  write.csv(res_mat, outfile, row.names = FALSE)
  log_message(paste("Results saved to:", outfile))
} else {
  log_message("No successful results to save")
}

# Optional: Check detailed logs
log_message(paste("Worker logs available in:", worker_logs_dir))