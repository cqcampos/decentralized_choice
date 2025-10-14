  ################################################################################
  # Script: likelihood_speed_test.R 
  # Author: Ryan Lee, Chris Campos  
  #
  # Loads in data for estimating the demand model 
  # Then, calculates likelihood and gradients for a chosen parameters
  # This script can be run in a Windows desktop, and meant for debugging likelihood code
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
  
  
  win_os <- 1
  if (win_os){
    dir <- "Z:/decentralized_choice"
    num_cores <- 1
  } else{
    # Set up number of cores to use (set 1 for windows PC)
    num_cores <- parallelly::availableCores() - 1
    dir <- "/Volumes/lausd/decentralized_choice"
  }

  setwd(dir)
  
  # Set this to path where the repo is located
  git_path <- file.path(Sys.getenv("HOME"), "Documents/GitHub/demand-estimation")

  # Log file setup
  log_file <- file.path(git_path, "estimation_log_runs.txt")
  
  
  # Load necessary libraries
  library(parallel)    # For parallel processing
  library(numDeriv)    # For numerical derivatives
  library(haven)       # For reading Stata files
  library(doParallel)  # For foreach parallel backend
  library(foreach)     # For parallel loops
  library(memoise)

  # Source helper files 
  source(file.path(git_path, "likelihood.R"))
  source(file.path(git_path, "get_data.R"))
  source(file.path(git_path, "objfn_no_partition.R"))
  source(file.path(git_path, "likelihood_mixture_types.R"))
  
  # Basic setup
  J <- 40                     # Number of schools
  R <- 3                     # Number of simulation draws for estimation
  R_post <- 1000              # Number of draws for calculating posteriors
  K <- 2                      # Number of mass points
  lambda <- 0.05              # Smoothing parameter
  homog <- 0                  # Treat charter schools as homogeneous (1) or heterogeneous (0)
  eta_out <- FALSE            # Is cost = exp() + eta or  exp(eta)?
  parallel_flag <- 0          # Use parallel processing? (1 = yes, 0 = no)
  Nblocks <- 1
  cl <- NULL
  
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
  start_time <- Sys.time()
  
  # We can parallelize the generation of simulation draws if needed
  sims_theta <- array(rnorm(N * R), dim = c(N, R))
  sims_c <- array(rnorm(N * J* R), dim = c(N, J, R))

  
  # Set initial parameter vector
  set.seed(123456)
  omega_choice_0 <- rnorm(n=42, mean=0,sd=0.02)
  omega_choice_0[1] <- -1.4
  omega_choice_0[14] <- -0.4
  omega_choice_0[28] <- 0.2
  omega_choice_0[29] <- -0.1
  omega_choice_0[28] <- 0.1
  
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
 
    omega_choice_0 <- c(
      -1.439757386, -1.924298342, -1.156092145, -1.383654371, -2.143312701,
      -1.261164888, -1.839703127, -2.107583895, -2.131663980, -0.934424440,
      -1.788556299, -1.699358249, -1.745399335, -0.366723408, -2.141872698,
      -1.313839011, -1.611662053, -1.107523843, -0.989861550, -1.958538161,
      -1.622796774, -1.258661378, -0.826065717, -1.775048373, -1.381252334,
      -1.986319901, -1.477205276, -2.152473038, -2.200929839,  0.199855041,
      -1.724229320, -1.802579615, -1.268705770, -1.047908769, -0.583168839,
      -1.846365769, -1.854806063, -1.181459273,  0.649342982, -1.293020795,
      -0.009803445, -0.245376369, -0.050943571, -0.324301425,  0.080543307,
      -0.216748835, -0.095833977,  0.143598405,  0.177957910, -0.022755289,
      0.551930736,  0.139519914,
      -0.232349101, -0.006922296,  0.032777485,  0.005049087, -0.066320875,
      0.004387379,  0.029886266,  0.023681391,  0.011126830, -0.007769884,
      -0.005079139,  0.079801487,  0.034247535,
      0.001411225, 
      -0.425049015, -0.002067020, -0.180884985,  0.181240487, -0.048737784,
      0.048204893,  0.008139680,  0.007193400, -0.065890145, -0.024599443,
      0.004106426, -0.298740757, -0.050776330,
      -0.388587959, # ln(sd(eta))
      -0.3, 0.3,   # ln(sd(theta_1)), ln(sd(theta_2))
       0, -0.5 #alpha_1, mu(theta_1)
    )
  }
  

  # omega_choice_0 <- omega_choice_0[80:84]
  
  # Full matrix to pass to the maximization routine:
  data_max <- list(K = K, Nblocks = Nblocks, A = A, A_dum = A_dum, Achoice=Achoice, Z = Z, E=E, X = X, D = D, W = W,
                   A_big = A_big, Z_big = Z_big, pi_i = pi_i, lambda = lambda, homog = homog,
                   choice_col = choice_col, E_col = E_col, eta_out = eta_out, num_cores = num_cores,
                   sims_theta = sims_theta, sims_c = sims_c, cl = cl, log_file = log_file)

  # source(file.path(git_path, "objfn_no_partition.R"))
  # Calculate analytic gradients and numerically approximated gradients
  # If the DGP in the likelihood calculation and the gradients calculation align,
  # these two types of gradients should be roughly the same
  #source(file.path(git_path, "likelihood_mixture_types.R"))
  
  #om1 <- c(-0.388587959, -0.3, 0.3, 0.232, -1)
  #om2 <- c(-0.388587959, -0.3, 0.3, (0.232+1e-3), -1)
  #q1 <- objfn_no_partition(om1, data_max, grad=FALSE)$Q
  #q2 <-objfn_no_partition(om2, data_max, grad=FALSE)$Q
  #print(q1)
  #print(q2)
  #(q2-q1)/1e-3
  

  
  # Some other data exploration  
  explore_data <- TRUE
  if (explore_data){
    
    source(file.path(git_path, "objfn_no_partition.R"))
  
    objfn_wrap <- function(par){
     objfn_no_partition(omega=par, data=data_max, grad=TRUE) 
    }
    
    mlmem <- memoise(objfn_wrap)
    
    fn <- function(par){
      mlmem(par)$Q
    }
    gr <- function(par){
      mlmem(par)$G
    }
    
    s2 <- Sys.time()
    
    # Maximization options (to be passed to optim)
    options_fmin <- list(maxit=1)
    
    optim_result <- optim(omega_choice_0,
                          fn = fn,
                          gr = gr,
                          method = "L-BFGS-B", 
                          lower  = -Inf,
                          upper  =  Inf,
                          control = options_fmin)
    
    
    
    print(Sys.time()-s2)
    
    # Calculate analytic gradients and numerically approximated gradients
    # If the DGP in the likelihood calculation and the gradients calculation align,
    s1 <- Sys.time()
    
    # Maximization options (to be passed to optim)
    optim_result <- optim(omega_choice_0,
                          fn = function(x) { objfn_no_partition(x, data_max, grad=FALSE)$Q },
                          gr = function(x) { objfn_no_partition(x, data_max, grad=TRUE)$G },
                          method = "L-BFGS-B", 
                          lower  = -Inf,
                          upper  =  Inf,
                          control = options_fmin)
    #num_grad_eta_in <- grad(function(om){objfn_no_partition(om, data_max, grad=FALSE)$Q}, omega_choice_0, method = "simple") 
    print(Sys.time()-s1)
    
    # these two types of gradients should be roughly the same
    #func_grad_eta_in <- objfn_no_partition(omega_choice_0, data_max, grad=TRUE)$G
    #func_grad_eta_in <- objfn_no_partition(omega_choice_0, data_max, grad=FALSE)$Q

    #num_grad_eta_in
    #func_grad_eta_in
    
    
    
    data_max$eta_out <- TRUE
 
    
    func_grad_eta_out <- objfn_no_partition(omega_choice_0, data_max, grad=TRUE)$G
    #num_grad_eta_out <- grad(function(om){objfn_no_partition(om, data_max, grad=FALSE)$Q}, omega_choice_0, method = "simple") 
    
    # num_grad_eta_out
    
    #func_grad_eta_out 
    
    
    
    if (FALSE){
    A_share <- colMeans(A)
    mean_utilities <- c(
      -1.566700359, -2.278722511, -1.463522005, -0.984785556, -2.149018592,
      -1.083422348, -1.341984054, -2.052743762, -2.194669979, -0.842057901,
      -1.797736056, -1.265049221, -1.818418941,  0.278898063, -2.189067360,
      -2.130892197, -1.686335032, -1.019400670, -0.935368860, -1.969482109,
      -1.488916811, -1.058695981, -0.603619599, -1.482008576, -1.249038789,
      -2.030492086, -1.762425437, -2.112614460, -2.291001626,  0.449211815,
      -1.803888168, -1.702881879, -1.321203466, -0.825357652, -0.761620354,
      -1.828773147, -1.913058623, -1.302095382,  0.142876902, -4.994268403
    )
    
    
    
    plot(y=A_share, x=mean_utilities, 
         xlab = "Estimated j-specific mean utility", 
         ylab = "Share of students in sample who applied to j")
  
    
    
    offer_accept <- matrix(0, nrow=J, ncol = 5)
    for (j in (1:J)){
      offer_j <- Z[,j]==1
      n_offer_j <- sum(offer_j)
      n_accept_j <- sum(offer_j == 1 & (E_col-1) == j)
      
      print(paste("Program:", j))
      print(paste("Num offer:", n_offer_j))
      print(paste("Num Accepted:", n_accept_j))
      print(paste("Share", n_accept_j/n_offer_j))
      
      
      offer_accept[j, 1] <- j
      offer_accept[j, 2] <- sum(choice_col == j)
      offer_accept[j, 3] <- n_offer_j
      offer_accept[j, 4] <- n_accept_j
      offer_accept[j, 5] <- n_accept_j/n_offer_j
      
    }
    
    df_share <- data.frame(offer_accept)
    names(df_share) <- c("j", "applied", "n_offer", "n_accept", "enroll_z_1")
    
    
    barplot(share~j, df_share, xlab = 1:J)
    }
  }