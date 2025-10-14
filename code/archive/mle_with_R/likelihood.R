################################################################################
# Script: likelihood.R 
# Author: Chris Campos, Ryan Lee, Connor Fogal
#
# Calculates likelihood contributions for demand model 
# 
# omega: vector of parameters
# data: matrix of data 
# grad: flag for calculating gradients
# worker_id: ID of the worker for logging (optional; for parallel processing)
# log_file: file path for logging output (optional)
# posterior: option to return posterior mean of theta for each individual 
################################################################################
likelihood <- function(omega, data, grad = FALSE, worker_id=0, log_file=NULL, posterior=FALSE) {
  # Set up logging
  if (is.null(log_file) && !is.null(data$log_file)) {
    log_file <- data$log_file
  }
  if (!is.null(data$worker_id)) {
    worker_id <- data$worker_id
  }
  # Worker-specific logging function
  lk_log <- function(msg) {
    if (!is.null(log_file)) {
      timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      formatted_msg <- paste0("[Worker ", worker_id, " (likelihood) - ", timestamp, "] ", msg)
      # We use cat with file argument to avoid locking issues
      cat(formatted_msg, "\n", file = log_file, append = TRUE)
    }
  }
  
  ##############################################################################
  ############################# Data Preparation ###############################
  ##############################################################################
  # Initialize outputs
  likelihood_i <- NULL
  G_i <- NULL
  
  
  # Extract data from list
  A <- data$A
  A_dum <- data$A_dum
  Achoice <- data$Achoice 
  Z <- data$Z
  X <- data$X
  E <- data$E
  W <- data$W
  A_big <- data$A_big
  Z_big <- data$Z_big
  D <- data$D
  pi_i <- data$pi_i
  lambda <- data$lambda
  homog <- data$homog
  sims_theta <- data$sims_theta
  sims_c <- data$sims_c
  
  # Determine dimensions
  Q_X <- ncol(X)
  Q_D <- dim(D)[2] 
  Q_W <- ncol(W)
  Q_A <- dim(A_big)[2] # Can apply to J + 1 portfolios 
  Q_Z <- dim(Z)[2] # Can receive offers from J programs
  N   <- dim(X)[1]
  J   <- dim(Z)[2] # J programs 
  R_draws <- dim(sims_theta)[2] 
  K <-  1 # Will update this later when moving to mixture model 

  ##############################################################################
  ############################ Extract Parameters ##############################
  ##############################################################################
  # First parameters correspond to school mean utilities (delta_j)
  l_param <- 1 
  # If only one choice mean utility, create J copies, else extract each mean utility 
  if(homog==1){
    delta_j <- rep(omega[l_param], J)
    l_param <- 2
  } else{
    delta_j <- omega[l_param:(l_param + J - 1)]
    l_param <- l_param + J
  }
  
  # Other parameters of the utility specification 
  beta_X <- matrix(omega[l_param:(l_param + Q_X - 1)], nrow = Q_X, ncol = 1)
  l_param <- l_param + Q_X
  psi_D <- matrix(omega[l_param:(l_param + Q_D - 1)], nrow = Q_D, ncol = 1)
  l_param <- l_param + Q_D
  
  # SD of eta_i (goes in costs)
  sigma_c <- exp(omega[l_param])
  l_param <- l_param + 1
  
  # SD of theta_i 
  sigma_theta <- exp(omega[l_param:(l_param + K - 1)])
  l_param <- l_param + K  
  
  # cost parameters (only need fixed cost)
  c_f <- omega[l_param:(l_param + Q_W - 1)]
  l_param <- l_param + Q_W
  
  # Parameters for the type distribution (this will change once we estimate types)
  mlogit_p_theta <- 1
  mu_theta <- 0
  
  if(posterior==TRUE){
    log_prior_r <- matrix(rep(0, N*R_draws), nrow=N, ncol=R_draws)
    theta_r <- matrix(rep(0, N*R_draws), nrow = N, ncol = R_draws)
  }
  
  ##############################################################################
  ########################### Calculations #####################################
  ##############################################################################
  
  ######################### Construct Mean Utilities ###########################
  # Construct mean utilities (U_bar) 
  # We have a constant delta_j (group_matrix) so let's make a matrix of that constant 
  # Mean utilities first 
  group_matrix <- matrix(rep(delta_j, each = N), nrow = N)  # N x J matrix
  # Next is the X*beta_X part 
  fixed_effect <- X %*% beta_X  # N x 1 vector
  # let's add them together and that's our U_bar (random coefficients come later)
  U_bar <- sweep(group_matrix, 1, fixed_effect, FUN = "+")
  # Add the distance (D) component (summing over Q_D for each school j)
  D_component <- matrix(0, nrow = N, ncol = J)

  # Takes the sum of each row of D[,,j]* psi_D which is just d_{ij}*X_i*psi_D
  # It is the multiplication by X_i that initially has D[,,,] 
  # After this, U_bar will be almost complete, just missing theta_i 
  for (j in 1:J) {
    D_component[, j] <- rowSums(D[,,j] * matrix(psi_D, nrow = N, ncol = Q_D, byrow = TRUE), na.rm=TRUE)
  }
  U_bar <- U_bar + D_component
  
  ############################ Application Costs ###############################
  applied <- apply(A, 1, sum)
  C_f_bar <- matrix(rep(exp(W %*% c_f), Q_A), ncol = Q_A)
  C_f_bar[,1] <- 0 

  ############################# Type Probabilities #############################
  #K==1 for now but will update with types 
  if (K == 1) {
    mlogit_p_i <- matrix(1, nrow = N)
    p_theta <- matrix(1, nrow = N)
    mlogit_p_0 <- 0
    p_theta_0 <- 1
  }
  
  # Useful information later
  hasOffers <- rowSums(Z) > 0
  applied <- rowSums(A) > 0
  enrolled <- rowSums(E[,2:Q_A])>0
  noEnrollmentStageChoice <- (applied ==0 | (applied==1 & hasOffers==0) )
  hasEnrollmentStageChoice <- (applied==1 & hasOffers==1)
  applied_chose_school <- (applied==1 & hasOffers==1 & enrolled==1)
  applied_chose_outside <- (applied==1 & hasOffers==1 & enrolled==0)
  
  ##############################################################################
  ############################# Simulated Likelihood ###########################
  ##############################################################################
  # Allocate arrays for simulated likelihood and its gradient
  Q_params <- length(omega)
  likelihood_sim <- array(0, dim = c(N, R_draws))
  dlikelihood_sim <- array(0, dim = c(N, Q_params, R_draws))
  
  # Progress tracking variables
  progress_interval <- max(1, floor(R_draws/10))  # Show progress every 10% or at least once
  start_time <- Sys.time()
  
  # Vectorized calculations for simulated likelihood
  for (r in 1:R_draws) {
    # Show progress
    if (r %% progress_interval == 0 || r == 1 || r == R_draws) {
      elapsed = (Sys.time() - start_time) 
      rate = elapsed / r 
      estimated_total = rate * R_draws
      estimated_remaining = estimated_total - elapsed
      
      lk_log(sprintf("Processing simulation draw %d/%d (%.1f%%) - Est. remaining: %.1f min", 
                     r, R_draws, 100*r/R_draws, as.numeric(estimated_remaining)/60))
    }
    
    ############################# Simulation Draws #############################
    # Extract simulation draws
    sims_theta_r <- matrix(sims_theta[, r], ncol = 1)
    sims_C_r <- sims_c[,r]  
    
    #This effectively imposes the constraint that the random coefficient is identical across all choice alternatives
    if (homog == 1) {
      sims_C_r <- matrix(rep(sims_C_r, J), ncol = J)
    }
    # Transform draws for theta and cost (eta)
    theta <- mu_theta + sigma_theta * sims_theta_r
    C_j <- sigma_c * sims_C_r  
    
    ###################### Finish cost calculations ############################
    eta_mat <- cbind( matrix(0, nrow=N, ncol=1), sigma_c * sims_C_r )  # N×(J+1)
    C_total <- C_f_bar * exp(eta_mat) # Random coefficient is in exponential 
    # make first column of costs zero (outside option)
    C_total <- cbind(matrix(0, nrow = N, ncol = 1), C_total[,2:ncol(C_total)])
    
    ############ Calculate app. stage expected utilities #######################
    # Create a matrix to keep track of Emax of different offer sets 
    # Dimensions: N × Q_Z, where Q_Z = J+1 (no offers + one per school)
    Emax <- matrix(0, nrow = N, ncol = J+1)
    U_tilde <- U_bar + matrix(rep(theta, J), nrow=N, ncol=J) 
    exp_U_0 <- exp(0)  # Outside option tau is just 0
    exp_U_1 <- exp(U_tilde)  # School options

    # Portfolio 1: Only outside option
    # Only option is outside option
    Emax[, 1] <-  0.577  # log(sum) + Euler's constant
    # All other portfolios: 
    for (j in 1:J) {
      # Options are: attend school j or choose outside option + Euler's constant
      Emax[, j+1] <- log( exp_U_0 + exp_U_1[, j] ) + 0.577
    }
    
    # With Emax can now calculate expected utilities (net of costs) for each offer set 
    e <- matrix(0, nrow = N, ncol = Q_A)
    # Portfolio 1: Not applying anywhere
    # Expected utility is just the utility of the outside option
    e[, 1] <- Emax[, 1]/lambda  # No applications -> only get the "no offers" outcome
    for (j in 1:J) {
      # Expected utility is weighted average of two possible outcomes:
      # 1. Get an offer (probability pi[,j]) -> utility is Emax[,j+1]
      # 2. Don't get an offer (probability 1-pi[,j]) -> utility is Emax[,1]
      # 3. and of course, subtract out application costs 
      e[, j+1] <- (pi_i[, j] * Emax[, j+1] + (1 - pi_i[, j]) * Emax[, 1] - C_total[, j+1])/lambda
    }
        # Exponentiate to then calculate probabilities 
    e <- exp(e)

    ######################## Application Probabilities #########################
    row_sums <- rowSums(e, na.rm=TRUE)
    P_A <- e / row_sums
    
    # Add outside option indicator 
    A_dum_full <- cbind(1-rowSums(A_dum), A_dum)    
    # Choose the correct P_A for observation i 
    P_A_i <- rowSums(P_A * A_dum_full, na.rm=TRUE)
    #P_A_i[P_A_i < 1e-6] <- 1e-6  # Rounding issues

    ######################### Enrollment Probabilities #########################
    # Only options with offers are considered
    numerator <- Z * exp_U_1
    # Create the denominator
    denominator <- 1 + rowSums(numerator)
    # Initialize probability matrix
    P_E <- matrix(0, nrow=N, ncol=J+1)
    # Identify students with no offers
    no_offers <- rowSums(Z) == 0
    # For students with no offers, outside option probability is 1
    P_E[no_offers, 1] <- 1
    # For students with offers, calculate according to the formula
    P_E[!no_offers, 1] <- 1 / denominator[!no_offers]  # Outside option
    for (j in 1:J) {
      P_E[!no_offers, j+1] <- numerator[!no_offers, j] / denominator[!no_offers]  # School j
    }
    # Probability of what student actually did 
    P_E_i <- rowSums(P_E * E, na.rm=TRUE)
    #P_E_i[P_E_i < 1e-6] <- 1e-6  # Rounding issues
  
    ################ Draw-specific likelihood contribution #####################
    likelihood_r <- P_A_i * P_E_i
    likelihood_sim[, r] <- likelihood_r
    
    # Two additional objects for posterior calculations
    if(posterior==TRUE){
      log_prior_r[,r]  <- -0.5 * (sims_theta_r^2)
      theta_r[,r] <- theta 
    }
    
    ##########################################################################
    ######################### Gradient Calculations ##########################
    ##########################################################################
    if (grad) {
      ########################## Initial Preparation ###########################
      # Initialize empty gradient matrices; we sum these at the end 
      dlog_P_A <- matrix(0, nrow = N, ncol = Q_params)
      dlog_P_E <- matrix(0, nrow = N, ncol = Q_params)

      q <- 1 # Index for parameters in the gradient
      
      # Create an array keeping track of all possible enrollment probabilities 
      # Probability choose j if face offer portfolio z (P(E_i=j|Z=z))
      P_E_Z <- array(0, dim = c(N, J+1, J+1))  # Dimensions: N individuals × Q_Z portfolios × J+1 options
      # For portfolio 1 (no offers), only outside option is available with probability 1 (all other probs are zero)
      P_E_Z[, 1, 1] <- 1
      # For each portfolio j (offer from school j only)
      for (j in 1:J) {
        # Calculate probability of choosing the outside option (option 1)
        P_E_Z[, j+1, 1] <- exp_U_0 / (exp_U_0 + exp_U_1[, j])
        
        # Calculate probability of choosing school j (option j+1)
        P_E_Z[, j+1, j+1] <- exp_U_1[, j] / (exp_U_0 + exp_U_1[, j])
      }
      # Keep track of application probabilities, excluding the outside option
      P_A_short <- P_A[, 2:(Q_A)]  # Exclude outside option (first column)
      
      #################### School mean utilities ##############################
      if(homog==1){
        dlog_P_A_temp <- matrix(0, nrow = N, ncol = 1)
        dlog_P_E_temp <- matrix(0, nrow = N, ncol = 1)
        ###### APPLICATION STAGE GRADIENT ###### 
        # # Calculate the dEmax for all potential offer portfolios and then assign correct one below
        dEmax <- matrix(0, nrow = N, ncol = J+1)
        # No effect on the "no offers" portfolio
        dEmax[, 1] <- 0
        for(j in 1:J) {
          # For portfolio j+1 (offer from school j), calculate derivative of log-sum (note pi_i missing but added below)
          dEmax[, j+1] <-  P_E_Z[, j+1, j+1]  
        }
       
        
        # Calculate de for application portfolios and dVE for enrollment gradient
        gA <- matrix(0, nrow = N, ncol = Q_A)
        gE <- matrix(0, nrow = N, ncol = Q_A)
        # No application portfolio not affected
        gA[, 1] <- 0
        gE[, 1] <- 0
        # For each school application
        for (a in 1:J) {
          # Expected derivative from offer outcomes
          gA[, a+1] <- pi_i[, a] * dEmax[, a+1] + (1 - pi_i[, a]) * dEmax[, 1]
          # Portion for the enrollment gradient
          gE[, a+1] <- dEmax[, a+1]
        }
        # Calculate gA
        # 1. pi*P_E for all applications (not including not applying)
        chosenA <- rowSums(gA[, 2:Q_A] * A_dum, na.rm=TRUE)
        # 2. sum_a P_A*pi*P_E 
        gA <- gA[, 2:Q_A]
        # 3. Gradient calculation 
        dlog_P_A_temp[applied] <- (chosenA[applied] - rowSums(P_A_short * gA,  na.rm=TRUE)[applied] )/lambda 
        dlog_P_A_temp[!applied] <- (0 - rowSums(P_A_short * gA,  na.rm=TRUE)[!applied] )/lambda 
        
        ###### ENROLLMENT STAGE GRADIENT ###### 
        # Calculate gE for enrollment stage
        chosenE <- rowSums(gE[, 2:Q_A] * A_dum, na.rm=TRUE)
        # Compute enrollment stage gradient by weighting deltaU by choice probabilities
        dlog_P_E_temp[applied_chose_school] <- (1- chosenE[applied_chose_school] ) # For those with offers and chose school
        dlog_P_E_temp[applied_chose_outside] <- -(chosenE[applied_chose_outside]) # For those with offers but chose outside option
        dlog_P_E_temp[noEnrollmentStageChoice] <- 0 # No info if do not have offers at this stage 
        
        # Assign to relevant column in the gradient matrix 
        dlog_P_A[,q] <- dlog_P_A_temp
        dlog_P_E[,q] <- dlog_P_E_temp
      }
      q <- q + 1

      #################### Covariate coefficients ##############################
      # Loop over each covariate (from 1 to Q_X)
      for(x in 1:Q_X) {
        
        ###### APPLICATION STAGE GRADIENT ###### 
        # Calculate dEmax directly for simplified portfolio structure
        dEmax <- matrix(0, nrow = N, ncol = J+1)
        Penroll <- matrix(0, nrow = N, ncol = J+1)
        # No effect on the "no offers" portfolio (outside option only)
        dEmax[, 1] <- 0
        # For each school offer portfolio j+1
        for (j in 1:J) {
          # For portfolio j+1 (offer from school j), calculate derivative of log-sum
          # Derivative equals probability of choosing school j times the covariate value (missing pi)
          # Note that this also works for the enrollment gradient 
          dEmax[, j+1] <- P_E_Z[, j+1, j+1] * X[, x]
          Penroll[, j+1] <- P_E_Z[, j+1, j+1]
        }
        # Calculate gA for application portfolios
        gA <- matrix(0, nrow = N, ncol = Q_A)
        # No application portfolio not affected 
        gA[, 1] <- 0
        
        # For each school application
        for (a in 1:J) {
          # Expected derivative from offer outcomes
          gA[, a+1] <- pi_i[, a] * dEmax[, a+1] + (1 - pi_i[, a]) * dEmax[, 1]
        }
        chosenA <- rowSums(gA[, 2:Q_A] * A_dum, na.rm=TRUE)
        # Gradient for application probability
        gA <- gA[, 2:Q_A]
        dlog_P_A[applied, q] <- (chosenA[applied] - rowSums(P_A_short * gA, na.rm=TRUE)[applied] )/lambda 
        dlog_P_A[!applied, q] <- (chosenA[!applied] - rowSums(P_A_short * gA, na.rm=TRUE)[!applied] )/lambda 
        
        ###### ENROLLMENT STAGE GRADIENT ###### 
        # Calculate gE for enrollment stage
        chosenE <- rowSums(Penroll[, 2:Q_A] * A_dum)
        # Compute enrollment stage gradient by weighting deltaU by choice probabilities
        dlog_P_E_temp[applied_chose_school] <- (1- chosenE[applied_chose_school])*X[applied_chose_school,x] # For those with offers and chose school
        dlog_P_E_temp[applied_chose_outside] <- -(chosenE[applied_chose_outside])*X[applied_chose_outside,x] # For those with offers but chose outside option
        dlog_P_E_temp[noEnrollmentStageChoice] <- 0 # No info if do not have offers at this stage 
        dlog_P_E[, q] <- dlog_P_E_temp
        q <- q + 1
      }

      
      #################### Distance coefficients ##############################
      for(d in 1:Q_D) {
        ###### APPLICATION STAGE GRADIENT ###### 
        dEmax <- matrix(0, nrow = N, ncol = J+1)
        Penroll <- matrix(0, nrow = N, ncol = J+1)
        # For each school offer portfolio j+1
        for (j in 1:J) {
          dEmax[, j+1] <- P_E_Z[, j+1, j+1] * D[, d, j]
          Penroll[, j+1] <- P_E_Z[, j+1, j+1]
        }
        
        # Calculate gA for application portfolios
        gA <- matrix(0, nrow = N, ncol = Q_A)
          
        # For each school application
        for (a in 1:J) {
          # Expected derivative from offer outcomes
          gA[, a+1] <- pi_i[, a] * dEmax[, a+1] + (1 - pi_i[, a]) * dEmax[, 1]
        }
        chosenA <- rowSums(gA[, 2:Q_A] * A_dum, na.rm=TRUE)
        # Gradient for application probability
        gA <- gA[, 2:Q_A]
        dlog_P_A[applied, q] <- (chosenA[applied] - rowSums(P_A_short * gA, na.rm=TRUE)[applied] )/lambda 
        dlog_P_A[!applied, q] <- (chosenA[!applied] - rowSums(P_A_short * gA, na.rm=TRUE)[!applied] )/lambda 
        ###### ENROLLMENT STAGE GRADIENT ###### 
        chosenE <- rowSums(Penroll[, 2:Q_A]*A_dum, na.rm=TRUE)
        relCovariate <- rowSums(D[,d,]*A_dum)
        dlog_P_E_temp[applied_chose_school] <- (1- chosenE[applied_chose_school])*relCovariate[applied_chose_school] # For those with offers and chose school
        dlog_P_E_temp[applied_chose_outside] <- -(chosenE[applied_chose_outside])*relCovariate[applied_chose_outside] # For those with offers but chose outside option
        dlog_P_E_temp[noEnrollmentStageChoice] <- 0 # No info if do not have offers at this stage 
        dlog_P_E[, q] <- dlog_P_E_temp
        q <- q + 1
      }
      
      ###################### SD of eta_i #####################################
      # Derivative of costs 
      dC <- C_total*sims_C_r[,1] # Note that this is N X (J +1)
      dC <- dC[, 2] # Keep one that corresponds to applying 
      dlog_P_A[applied, q] <- -(dC[applied] - rowSums(P_A_short * dC, na.rm=TRUE)[applied])/lambda 
      dlog_P_A[!applied, q] <- (rowSums(P_A_short * dC, na.rm=TRUE)[!applied] )/lambda 
      q <- q + 1
      ###################### SD of theta_i ####################################
      for (k in 1:K) {
        ###### APPLICATION STAGE GRADIENT ###### 
        dU <- sims_theta_r 
        # Calculate dEmax directly for simplified portfolio structure
        dEmax <- matrix(0, nrow = N, ncol = J+1)
        Penroll <- matrix(0, nrow = N, ncol = J+1)
        # No effect on the "no offers" portfolio (outside option only)
        dEmax[, 1] <- 0
        # For each school offer portfolio j+1
        for (j in 1:J) {
          # For portfolio j+1 (offer from school j), calculate derivative of log-sum
          # Derivative equals probability of choosing school j times the covariate value (missing pi)
          # Note that this also works for the enrollment gradient 
          dEmax[, j+1] <- P_E_Z[, j+1, j+1] * dU
          Penroll[, j+1] <- P_E_Z[, j+1, j+1]
        }
        # Calculate gA for application portfolios
        gA <- matrix(0, nrow = N, ncol = Q_A)
        # No application portfolio not affected 
        gA[, 1] <- 0
        # For each school application
        for (a in 1:J) {
          # Expected derivative from offer outcomes
          gA[, a+1] <- pi_i[, a] * dEmax[, a+1] + (1 - pi_i[, a]) * dEmax[, 1]
        }
        chosenA <- rowSums(gA[, 2:Q_A] * A_dum, na.rm=TRUE)
        # Gradient for application probability
        gA <- gA[, 2:Q_A]
        dlog_P_A[applied, q] <- (chosenA[applied] - rowSums(P_A_short * gA, na.rm=TRUE)[applied] )/lambda 
        dlog_P_A[!applied, q] <- (chosenA[!applied] - rowSums(P_A_short * gA, na.rm=TRUE)[!applied] )/lambda 
        ###### ENROLLMENT STAGE GRADIENT ###### 
        # Calculate gE for enrollment stage
        chosenE <- rowSums(Penroll[, 2:Q_A] * A_dum, na.rm=TRUE)
        # Compute enrollment stage gradient by weighting deltaU by choice probabilities
        dlog_P_E_temp[applied_chose_school] <- (1- chosenE[applied_chose_school])*dU[applied_chose_school] # For those with offers and chose school
        dlog_P_E_temp[applied_chose_outside] <- -(chosenE[applied_chose_outside])*dU[applied_chose_outside] # For those with offers but chose outside option
        dlog_P_E_temp[!hasEnrollmentStageChoice] <- 0 # No info if do not have offers at this stage 
        dlog_P_E[, q] <- dlog_P_E_temp
        q <- q + 1
      }

 
      ################## Application fixed cost coefficients ###################
      for (w in 1:Q_W) {
        # Derivative of costs 
        dC <- C_total*W[,w] # Note that this is N X (J +1)
        dC <- dC[, 2] # Keep one that corresponds to applying 
        dlog_P_A[applied, q] <- -(dC[applied] - rowSums(P_A_short * dC, na.rm=TRUE)[applied])/lambda 
        dlog_P_A[!applied, q] <- (rowSums(P_A_short * dC, na.rm=TRUE)[!applied] )/lambda 
        q <- q + 1
      }

      ##################### Compute gradiet for draw r #########################
      # Compute total gradient for this simulation draw
      dlikelihood_sim[,,r] <- (dlog_P_A + dlog_P_E) * likelihood_r  
      # Log an update on completion if at certain percentage 
      if (r %% progress_interval == 0 || r == R_draws) {
        lk_log(sprintf("Completed gradient calculations for draw %d/%d", r, R_draws))
      }
    }
  }
  # Final progress report
  total_time <- Sys.time() - start_time
  lk_log(sprintf("Completed all %d simulation draws in %.2f minutes", 
                 R_draws, as.numeric(total_time)/60))

  # Calculate final likelihood by averaging over simulation draws
  likelihood_i <- rowMeans(likelihood_sim, na.rm=TRUE)
  print(summary(likelihood_i))
  print(quantile(likelihood_i, probs=c(0.01, 0.05, 0.10, 0.15, 0.20, 0.40, 0.5, 0.60, 0.80, 0.90, 0.95, 0.99), na.rm=TRUE))

  # Ensure no zeros (could cause log(0) issues
  threshold_count <- sum(likelihood_i < 1e-6, na.rm=TRUE)
  #likelihood_i[likelihood_i < 1e-16] <- 1e-16

  percent_threshold <- round(100*threshold_count/length(likelihood_i), 2)
  lk_log(sprintf("Low likelihood warning: %d observations (%g%%) hit minimum threshold (1e-6)", 
                 threshold_count, percent_threshold))
  lk_log(sprintf("Pre-calculation check: %d NaNs in likelihood_i", sum(is.na(likelihood_i))))
  lk_log(sprintf("Likelihood summary - Min: %.3e, Q1: %.3e, Median: %.3e, Mean: %.3e, Q3: %.3e, Max: %.3e",
                 min(likelihood_i, na.rm=TRUE), quantile(likelihood_i, 0.25, na.rm=TRUE), median(likelihood_i, na.rm=TRUE),
                 mean(likelihood_i, na.rm=TRUE), quantile(likelihood_i, 0.75, na.rm=TRUE), max(likelihood_i, na.rm=TRUE)))
  
  # Calculate final gradient if requested
  if (grad) {
    # Average gradient over simulation draws
    dlikelihood_i <- rowMeans(dlikelihood_sim, dims = 2, na.rm = TRUE)
    print(summary(dlikelihood_i))
    
    # Calculate gradient of log likelihood
    likelihood_i_mat <- matrix(likelihood_i, nrow = length(likelihood_i) , ncol = ncol(dlikelihood_i), byrow = FALSE)
    G_i <- dlikelihood_i / likelihood_i_mat # N x No. Parameters
    #print(quantile(G_i, probs=c(0.05, 0.10, 0.20, 0.40, 0.5, 0.60, 0.80, 0.90, 0.95), na.rm=TRUE))
    
    # Check for NaN or Inf values
    count_gi_nan <- sum(is.nan(G_i), na.rm = TRUE)
    lk_log(sprintf("Pre-calculation check: %d NaNs in gradient", count_gi_nan  ))
  }
  lk_log("Likelihood function completed")
  
  theta_posterior_mu <- NULL
  theta_posterior_sd <- NULL
  if(posterior==TRUE){
    # Calculate posterior weights
    log_lik <- log(likelihood_sim)
    log_weights <- log_lik + log_prior_r
    # Normalize weights
    weight <- exp(log_weights)
    weights <- weight / rowSums(weight)
    
    # Calculate posterior means and variances
    theta_posterior_mu <- rowSums(weights * theta_r)
    theta_posterior_var <- rowSums(weights * (theta_r - theta_posterior_mu)^2)
    theta_posterior_sd <- sqrt(theta_posterior_var)
    
    
  }
  return(list(L = likelihood_i, G = G_i, 
              posterior_mu = theta_posterior_mu, 
              posterior_sd = theta_posterior_sd))
}




likelihood_fast_par <- function(omega, data, grad = FALSE, log_file=NULL, posterior=FALSE) {
  # Set up logging
  
  if (is.null(log_file) && !is.null(data$log_file)) {
    log_file <- data$log_file
  }
  
  
  # Worker-specific logging function
  lk_log <- function(msg) {
    if (!is.null(log_file)) {
      timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      formatted_msg <- paste0(" (likelihood) - ", timestamp, "] ", msg)
      # We use cat with file argument to avoid locking issues
      cat(formatted_msg, "\n", file = log_file, append = TRUE)
    }
  }
  
  ##############################################################################
  ############################# Data Preparation ###############################
  ##############################################################################
  # Initialize outputs
  likelihood_i <- NULL
  G_i <- NULL

  
  # Extract data from list
  A <- data$A
  A_dum <- data$A_dum
  Achoice <- data$Achoice 
  applied <- data$applied
  Z <- data$Z
  X <- data$X
  E <- data$E
  W <- data$W
  A_big <- data$A_big
  Z_big <- data$Z_big
  pi_i <- data$pi_i
  lambda <- data$lambda
  homog <- data$homog
  sims_theta <- data$sims_theta
  sims_c <- data$sims_c
  cl <- data$cl
  choice_col <- data$choice_col
  E_col <- data$E_col
  eta_out <- data$eta_out
  n_cores <- data$num_cores
  rm(data)
  gc()
  
  # Determine dimensions
  Q_X <- ncol(X)
  Q_D <- dim(D)[2] 
  Q_W <- ncol(W)
  Q_A <- dim(A_big)[2] # Can apply to J + 1 portfolios 
  Q_Z <- dim(Z)[2] # Can receive offers from J programs
  N   <- dim(X)[1]
  J   <- dim(Z)[2] # J programs 
  R_draws <- dim(sims_theta)[2] 
  K <-  1 # Will update this later when moving to mixture model 
  
  ##############################################################################
  ############################ Extract Parameters ##############################
  ##############################################################################
  # First parameters correspond to school mean utilities (delta_j)
  l_param <- 1 
  # If only one choice mean utility, create J copies, else extract each mean utility 
  if(homog==1){
    delta_j <- rep(omega[l_param], J)
    l_param <- 2
  } else{
    delta_j <- omega[l_param:(l_param + J - 1)]
    l_param <- l_param + J
  }
  
  # Other parameters of the utility specification 
  beta_X <- matrix(omega[l_param:(l_param + Q_X - 1)], nrow = Q_X, ncol = 1)
  l_param <- l_param + Q_X
  psi_D <- matrix(omega[l_param:(l_param + Q_D - 1)], nrow = Q_D, ncol = 1)
  l_param <- l_param + Q_D
  
  # SD of eta_i (goes in costs)
  sigma_c <- exp(omega[l_param])
  l_param <- l_param + 1
  
  # SD of theta_i 
  sigma_theta <- exp(omega[l_param:(l_param + K - 1)])
  l_param <- l_param + K  
  
  # cost parameters (only need fixed cost)
  c_f <- omega[l_param:(l_param + Q_W - 1)]
  l_param <- l_param + Q_W
 
  
  
  # Parameters for the type distribution (this will change once we estimate types)
  mlogit_p_theta <- 1
  mu_theta <- 0
  
  
  log_prior_r <- matrix(0, nrow=N, ncol=R_draws)
  theta_r <- matrix(0, nrow = N, ncol = R_draws)
  
  
  ##############################################################################
  ########################### Calculations #####################################
  ##############################################################################
  
  ######################### Construct Mean Utilities ###########################
  # Construct mean utilities (U_bar) 
  # We have a constant delta_j (group_matrix) so let's make a matrix of that constant 
  # Mean utilities first 
  group_matrix <- matrix(rep(delta_j, each = N), nrow = N)  # N x J matrix
  # Next is the X*beta_X part 
  fixed_effect <- X %*% beta_X  # N x 1 vector
  # let's add them together and that's our U_bar (random coefficients come later)
  U_bar <- sweep(group_matrix, 1, fixed_effect, FUN = "+")
  # Add the distance (D) component (summing over Q_D for each school j)
  D_component <- matrix(0, nrow = N, ncol = J)
  
  # Takes the sum of each row of D[,,j]* psi_D which is just d_{ij}*X_i*psi_D
  # It is the multiplication by X_i that initially has D[,,,] 
  # After this, U_bar will be almost complete, just missing theta_i 
  for (j in 1:J) {
    D_component[, j] <- rowSums(D[,,j] * matrix(psi_D, nrow = N, ncol = Q_D, byrow = TRUE), na.rm=TRUE)
  }
  U_bar <- U_bar + D_component
  
  ############################ Application Costs ###############################
  applied <- apply(A, 1, sum)
  C_f_bar <- matrix(rep(exp(W %*% c_f), Q_A), ncol = Q_A)
  C_f_bar[,1] <- 0 
  
  ############################# Type Probabilities #############################
  #K==1 for now but will update with types 
  if (K == 1) {
    mlogit_p_i <- matrix(1, nrow = N)
    p_theta <- matrix(1, nrow = N)
    mlogit_p_0 <- 0
    p_theta_0 <- 1
  }
  
  # Useful information later
  hasOffers <- rowSums(Z) > 0
  applied <- rowSums(A) > 0
  enrolled <- rowSums(E[,2:Q_A])>0
  noEnrollmentStageChoice <- (applied ==0 | (applied==1 & hasOffers==0) )
  hasEnrollmentStageChoice <- (applied==1 & hasOffers==1)
  applied_chose_school <- (applied==1 & hasOffers==1 & enrolled==1)
  applied_chose_outside <- (applied==1 & hasOffers==1 & enrolled==0)
  
  ##############################################################################
  ############################# Simulated Likelihood ###########################
  ##############################################################################
  # Allocate arrays for simulated likelihood and its gradient
  Q_params <- length(omega)
  likelihood_sim <- array(0, dim = c(N, R_draws))
  dlikelihood_sim <- array(0, dim = c(N, Q_params, R_draws))
  
  # Progress tracking variables
  progress_interval <- max(1, floor(R_draws/10))  # Show progress every 10% or at least once
  start_time <- Sys.time()
  
  # Helper function for calculating likelihood for simulation draw `r`
  l_and_grad_r <- function(r){
    
    
    ############################# Simulation Draws #############################
    # Extract simulation draws
    sims_theta_r <- matrix(sims_theta[, r], ncol = 1)
    # sims_C_r <- sims_c[,,r]  

    
    # Transform draws for theta and cost (eta)
    theta <- mu_theta + sigma_theta * sims_theta_r
    cur_eta_ij2 <- matrix(sims_c[,,r]^2, nrow=N, ncol=J)
    eta_ij <- sigma_c * sims_c[,,r]  
    
    ###################### Finish cost calculations ############################
    #eta_mat <- cbind( matrix(0, nrow=N, ncol=1), sigma_c * sims_C_r )  # N×(J+1)
    eta_mat <- cbind(matrix(0, nrow=N, ncol=1), eta_ij )
    
    if (eta_out){
      C_total <- C_f_bar + eta_mat
    } else{
      C_total <- C_f_bar * exp(eta_mat) # Random coefficient is in exponential 
    }
    rm(eta_mat)
    gc()
    
    # make first column of costs zero (outside option)
    C_total <- cbind(matrix(0, nrow = N, ncol = 1), C_total[,2:ncol(C_total)])
    
    ############ Calculate app. stage expected utilities #######################
    # Create a matrix to keep track of Emax of different offer sets 
    # Dimensions: N × Q_Z, where Q_Z = J+1 (no offers + one per school)
    Emax <- matrix(0, nrow = N, ncol = J+1)
    U_tilde <- U_bar + matrix(rep(theta, J), nrow=N, ncol=J) 
    exp_U_0 <- exp(0)  # Outside option tau is just 0
    exp_U_1 <- exp(U_tilde)  # School options
    
    # Portfolio 1: Only outside option
    # Only option is outside option
    Emax[, 1] <-  0.5772156649  # log(sum) + Euler's constant
    # All other portfolios: 
    for (j in 1:J) {
      # Options are: attend school j or choose outside option + Euler's constant
      Emax[, j+1] <- log( exp_U_0 + exp_U_1[, j] ) + 0.5772156649
    }
    
    
    # With Emax can now calculate expected utilities (net of costs) for each offer set 
    e <- matrix(0, nrow = N, ncol = Q_A)
    # Portfolio 1: Not applying anywhere
    # Expected utility is just the utility of the outside option
    e[, 1] <- Emax[, 1]/lambda  # No applications -> only get the "no offers" outcome
    for (j in 1:J) {
      # Expected utility is weighted average of two possible outcomes:
      # 1. Get an offer (probability pi[,j]) -> utility is Emax[,j+1]
      # 2. Don't get an offer (probability 1-pi[,j]) -> utility is Emax[,1]
      # 3. and of course, subtract out application costs 
      e[, j+1] <- (1/lambda)* pi_i[, j] * (Emax[, j+1])  + 
                  (1/lambda)*(1-pi_i[, j]) *(Emax[, 1]) - 
                  (C_total[, j+1])/lambda
    }
    # Exponentiate to then calculate probabilities 
    e <- exp(e)
    
    ######################## Application Probabilities #########################
    row_sums <- rowSums(e, na.rm=TRUE)
    P_A <- e / row_sums
    
    # Add outside option indicator 
    A_dum_full <- cbind(1-rowSums(A_dum), A_dum)    
    # Choose the correct P_A for observation i 
    P_A_i <- rowSums(P_A * A_dum_full, na.rm=TRUE)
    #P_A_i[P_A_i < 1e-6] <- 1e-6  # Rounding issues
    
    ######################### Enrollment Probabilities #########################
    # Only options with offers are considered
    numerator <- Z * exp_U_1
    # Create the denominator
    denominator <- 1 + rowSums(numerator)
    # Initialize probability matrix
    P_E <- matrix(0, nrow=N, ncol=J+1)
    # Identify students with no offers
    no_offers <- rowSums(Z) == 0
    # For students with no offers, outside option probability is 1
    P_E[no_offers, 1] <- 1
    # For students with offers, calculate according to the formula
    P_E[!no_offers, 1] <- 1 / denominator[!no_offers]  # Outside option
    for (j in 1:J) {
      P_E[!no_offers, j+1] <- numerator[!no_offers, j] / denominator[!no_offers]  # School j
    }
    # Probability of what student actually did 
    P_E_i <- rowSums(P_E * E, na.rm=TRUE)
    #P_E_i[P_E_i < 1e-6] <- 1e-6  # Rounding issues
    
    ################ Draw-specific likelihood contribution #####################
    likelihood_r <- P_A_i * P_E_i
  
    # Two additional objects for posterior calculations
    if(posterior==TRUE){
      #print("Printing log prior dimension")
      #print(dim(log_prior_r))
      #print("Printing sims_theta_r dimension")
      #print(dim(sims_theta_r))
      log_prior_r[,r]  <- (-0.5) * (sims_theta_r^2) 
      theta_r[,r] <- theta 
    }
    
    ##########################################################################
    ######################### Gradient Calculations ##########################
    ##########################################################################
    if (grad) {
      ########################## Initial Preparation ###########################
      
      # Create an array keeping track of all possible enrollment probabilities 
      # Probability choose j if face offer portfolio z (P(E_i=j|Z=z))
      P_E_Z <- array(0, dim = c(N, J+1, 2))  # Dimensions: N individuals × Q_Z portfolios × J+1 options
      # For portfolio 1 (no offers), only outside option is available with probability 1 (all other probs are zero)
      P_E_Z[, 1, 1] <- 1
      # For each portfolio j (offer from school j only)
      for (j in 1:J) {
        # Calculate probability of choosing the outside option (option 1)
        P_E_Z[, j+1, 1] <- exp_U_0 / (exp_U_0 + exp_U_1[, j])
        
        # Calculate probability of choosing school j (option j+1)
        P_E_Z[, j+1, 2] <- exp_U_1[, j] / (exp_U_0 + exp_U_1[, j])
      }
      
      # Keep track of application probabilities, excluding the outside option
      P_A_short <- P_A[, 2:(Q_A)]  # Exclude outside option (first column)
      
      
      ################# Application and Enrollment gradients ###################
      
      # Calculate de for application portfolios and dVE for enrollment gradient
      gA <- matrix(0, nrow = N, ncol = Q_A)
      gE <- matrix(0, nrow = N, ncol = Q_A)
      # No application portfolio not affected
      gA[, 1] <- 0
      gE[, 1] <- 0
      # For each application, calculate:
      for (a in 1:J) {
        # Derivative of E(u(max) \in offer)) w.r.t \delta_j
        # But it's part of the application gradients w.r.t other parameters
        gA[, a+1] <- pi_i[, a] * P_E_Z[, a+1, 2] 
        
        # Portion for the enrollment gradient
        gE[, a+1] <- P_E_Z[, a+1, 2]
      }
      
      # Calculate gA
      # 1. pi*P_E for all applications (not including not applying)
      chosenA <- rowSums(gA[, 2:Q_A] * A_dum, na.rm=TRUE)
      # 2. sum_a P_A*pi*P_E 
      gA <- gA[, 2:Q_A]
      # 3. Gradient calculation 
      
      # Calculate gE for enrollment stage
      chosenE <- rowSums(gE[, 2:Q_A] * A_dum, na.rm=TRUE)
      
      
      if(homog==1){
        
        # Application stage
        dlog_P_A_delta <- (chosenA - rowSums(P_A_short * gA,  na.rm=TRUE))/lambda
        
        # Compute enrollment stage gradient by weighting deltaU by choice probabilities
        #dlog_P_E_temp[applied_chose_school] <- (1- chosenE[applied_chose_school] ) # For those with offers and chose school
        #dlog_P_E_temp[applied_chose_outside] <- -(chosenE[applied_chose_outside]) # For those with offers but chose outside option
        #dlog_P_E_temp[noEnrollmentStageChoice] <- 0 # No info if do not have offers at this stage
        dlog_P_E_delta <- (ifelse(applied_chose_school, 1, 0) - chosenE)*(1-noEnrollmentStageChoice)
      } else{
        dlog_P_A_delta <- matrix(0, nrow = N, ncol = J)
        dlog_P_E_delta <- matrix(0, nrow = N, ncol = J)
        for (j in 1:J){
          # Application stage
          # ChosenA = pi_{i,a_i}*{P_E,a_i}
          dlog_P_A_delta[,j] <- ((1/lambda)*chosenA*(choice_col==j)) - ((1/lambda)*P_A_short[,j] * gA[, j])
          
          # Compute enrollment stage gradient by weighting deltaU by choice probabilities
          dlog_P_E_delta[,j] <- (applied_chose_school - chosenE)*(choice_col == j)*(1-noEnrollmentStageChoice)
        }
      } 
      
      #################### Covariate coefficients ##############################
      # Loop over each covariate (from 1 to Q_X)
      dlog_P_A_X <- matrix(0, nrow = N, ncol = Q_X)
      dlog_P_E_X <- matrix(0, nrow = N, ncol = Q_X) 
      for(x in 1:Q_X) {
        # App stage grad w.r.t x
        dlog_P_A_X[, x] <- ((1/lambda)*chosenA*X[,x]) - rowSums( (1/lambda)*gA * X[,x]*P_A_short ,  na.rm=TRUE)
        
        # enrollment stage grad
        dlog_P_E_X[, x] <- (ifelse(applied_chose_school, 1, 0) - chosenE)*X[,x]*(1-noEnrollmentStageChoice)
      }
      
      
      #################### Distance coefficients ##############################
      dlog_P_A_D <- matrix(0, nrow = N, ncol = Q_D)
      dlog_P_E_D <- matrix(0, nrow = N, ncol = Q_D) 
      
      # Loop over interaction between covariates and distances 
      for(d in 1:Q_D) {
        relCovariate <- rowSums(D[,d,]*A_dum)
        
        # App stage grad w.r.t x
        dlog_P_A_D[, d] <- ((1/lambda)*chosenA*relCovariate)- (rowSums((1/lambda)*gA * D[,d,]*P_A_short ,  na.rm=TRUE))
        
        # enrollment stage grad
        dlog_P_E_D[, d] <- (ifelse(applied_chose_school, 1, 0) - chosenE)*relCovariate*(1-noEnrollmentStageChoice)
      }
      
      ###################### SD of eta_i #####################################
      # Derivative of costs 
      #dC <- C_total*sims_C_r[,1] # Note that this is N X (J +1)
      #dC <- dC[, 2] # Keep one that corresponds to applying 
      #dlog_P_A_sd_eta <- -(ifelse(applied, dC, 0) - (rowSums(P_A_short * dC, na.rm=TRUE)))/lambda
      
      
      if (eta_out){
        # New cost gradient with heterogeneous eta_ij 
        dC <- eta_ij 
      } else{
        # When eta is in exp() 
        dC <- C_total[,2:ncol(C_total)] * eta_ij # Note that this is N X (J +1)
      }
      rm(eta_ij)
      gc()
      chosendC <- rowSums(dC * A_dum, na.rm=TRUE)
      dlog_P_A_sd_eta <- -((1/lambda)*chosendC - (rowSums((1/lambda)*P_A_short * dC, na.rm=TRUE)))

      
      ###################### SD of theta_i ####################################
      for (k in 1:K) {
        ###### APPLICATION STAGE GRADIENT ###### 
        dU <- theta
        dlog_P_A_theta <- ((1/lambda)*chosenA*dU -  rowSums( (1/lambda) * P_A_short * gA * dU[, rep(1, ncol(gA))], na.rm=TRUE))
        dlog_P_E_theta <- (ifelse(applied_chose_school, dU, 0) - chosenE*dU)*(1-noEnrollmentStageChoice)
      }
      
      ################## Application fixed cost coefficients ###################
      dlog_P_A_W <- matrix(0, nrow = N, ncol = Q_W)
      
      
      # This code can be shorter by putting if (eta_out) inside the w loop, but didn't want
      # the code to check eta_out==TRUE Q_W times
      if (eta_out){        
        for (w in 1:Q_W) {
          # Derivative of costs 
          dC <- C_f_bar[,2:Q_A]*W[,w] # Note that this is N X (J +1)
          dCapplied <- rowSums(dC*A_dum, na.rm=TRUE) # Keep one that corresponds to applied
          dlog_P_A_W[, w] <- -( (1/lambda)*dCapplied - rowSums((1/lambda)*P_A_short * dC, na.rm=TRUE))
        }
      } else{
        for (w in 1:Q_W) {
          # Derivative of costs 
          dC <- C_total*W[,w] # Note that this is N X (J +1)
          dC <- dC[, -1] # Keep ones that corresponds to applying 
          
          dCapplied <- rowSums(dC*A_dum, na.rm=TRUE) # Keep one that corresponds to applied
          dlog_P_A_W[, w] <- -( (1/lambda)*dCapplied - rowSums((1/lambda)*P_A_short * dC, na.rm=TRUE))
        }
      }
      
      # Matrix of 0 for matching dimension for matrix addition
      dlog_P_E_sd_eta <- matrix(0, nrow = N, ncol =1)
      
      dlog_P_E_W <- matrix(0, nrow = N, ncol = Q_W) 
      
      ## Combine all gradients, in order they should be combined
      dlog_P_A <- cbind(dlog_P_A_delta, dlog_P_A_X, dlog_P_A_D, dlog_P_A_sd_eta, dlog_P_A_theta, dlog_P_A_W)
      
      dlog_P_E <- cbind(dlog_P_E_delta, dlog_P_E_X, dlog_P_E_D, dlog_P_E_sd_eta, dlog_P_E_theta, dlog_P_E_W)
      
      
      ##################### Compute gradient for draw r #########################
      # Compute total gradient for this simulation draw
      dlikelihood_sim_r <- (dlog_P_A + dlog_P_E) * likelihood_r  
    }
    
    if (grad){
      return(list(likelihood_r=likelihood_r, dlikelihood_sim_r = dlikelihood_sim_r))
    } else {
      return(list(likelihood_r=likelihood_r, log_prior_r = log_prior_r, theta_r=theta_r ))
    }
    
  }
  
  res <- mclapply(1:R_draws, function(i){l_and_grad_r(i)}, mc.cores=n_cores)


  
  # Create the likelihood matrix (cbind of 1st list elements from each sublist)
  
  if (grad){
    likelihood_sim <- do.call(cbind, lapply(res, `[[`, 1))
    
    # Create the 3D gradient array from the 2nd list elements
    array_dim <- c(dim(res[[1]][[2]]), length(res))  # N, Q, 30
    dlikelihood_sim <- array(unlist(lapply(res, `[[`, 2)), dim = array_dim)
  } else{
    likelihood_sim <- do.call(cbind, lapply(res, `[[`, 1))
    log_prior_r <- do.call(cbind, lapply(res, `[[`, 2))
    theta_r <- do.call(cbind, lapply(res, `[[`, 3))
  }
  
  
  # Final progress report
  total_time <- Sys.time() - start_time
  lk_log(sprintf("Completed all %d simulation draws in %.2f minutes", 
                 R_draws, as.numeric(total_time)/60))
  
  # Calculate final likelihood by averaging over simulation draws
  #print(dim(likelihood_sim))
  # print(length(likelihood_sim))
  #print(likelihood_sim)
  print(dim(likelihood_sim))
  print(summary(likelihood_sim))
  likelihood_i <- rowMeans(likelihood_sim, na.rm=TRUE)
  print(summary(likelihood_i))
  print(quantile(likelihood_i, probs=c(0.01, 0.05, 0.10, 0.15, 0.20, 0.40, 0.5, 0.60, 0.80, 0.90, 0.95, 0.99), na.rm=TRUE))
  
  # Ensure no zeros (could cause log(0) issues
  threshold_count <- sum(likelihood_i < 1e-6, na.rm=TRUE)
  #likelihood_i[likelihood_i < 1e-16] <- 1e-16
  
  percent_threshold <- round(100*threshold_count/length(likelihood_i), 2)
  lk_log(sprintf("Low likelihood warning: %d observations (%g%%) hit minimum threshold (1e-6)", 
                 threshold_count, percent_threshold))
  lk_log(sprintf("Pre-calculation check: %d NaNs in likelihood_i", sum(is.na(likelihood_i))))
  lk_log(sprintf("Likelihood summary - Min: %.3e, Q1: %.3e, Median: %.3e, Mean: %.3e, Q3: %.3e, Max: %.3e",
                 min(likelihood_i, na.rm=TRUE), quantile(likelihood_i, 0.25, na.rm=TRUE), median(likelihood_i, na.rm=TRUE),
                 mean(likelihood_i, na.rm=TRUE), quantile(likelihood_i, 0.75, na.rm=TRUE), max(likelihood_i, na.rm=TRUE)))
  
  # Calculate final gradient if requested
  if (grad) {
    # Average gradient over simulation draws
    dlikelihood_i <- rowMeans(dlikelihood_sim, dims = 2, na.rm = TRUE)
    print(summary(dlikelihood_i))
    
    # Calculate gradient of log likelihood
    likelihood_i_mat <- matrix(likelihood_i, nrow = length(likelihood_i) , ncol = ncol(dlikelihood_i), byrow = FALSE)
    G_i <- dlikelihood_i / likelihood_i_mat # N x No. Parameters
    #print(quantile(G_i, probs=c(0.05, 0.10, 0.20, 0.40, 0.5, 0.60, 0.80, 0.90, 0.95), na.rm=TRUE))
    
    # Check for NaN or Inf values
    count_gi_nan <- sum(is.nan(G_i), na.rm = TRUE)
    lk_log(sprintf("Pre-calculation check: %d NaNs in gradient", count_gi_nan  ))
  }
  lk_log("Likelihood function completed")
  
  theta_posterior_mu <- NULL
  theta_posterior_sd <- NULL
  
  if(posterior==TRUE){
    # Calculate posterior weights
    print("Printing log_prior_r")
    print(log_prior_r)
    log_lik <- log(likelihood_sim)
    print("Printing log_lik")
    print(summary(log_lik))
    log_weights <- log_lik + log_prior_r
    
    
    print("Printing log_weights")
    print(summary(log_weights))
    # Normalize weights
    
    weights <- exp(log_weights)
    weights <- weights / rowSums(weights)
    print("Printing soft max weights")
    print(summary(weights))
    print("Printing theta_r")
    print(summary(theta_r))
    # Calculate posterior means and variances
    theta_posterior_mu <- rowSums(weights * theta_r)
    theta_posterior_var <- rowSums(weights * (theta_r - theta_posterior_mu)^2)
    theta_posterior_sd <- sqrt(theta_posterior_var)
    
    
  }
  return(list(L = likelihood_i, G = G_i, 
              posterior_mu = theta_posterior_mu, 
              posterior_sd = theta_posterior_sd))
}

