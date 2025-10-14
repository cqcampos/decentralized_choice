likelihood_K_mixture <- function(omega, data, grad = FALSE, log_file=NULL, posterior=FALSE) {
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
  
  if (FALSE){
    data <- data_max
    omega <- omega_choice_0 
    grad <- TRUE
    posterior <- FALSE
    homog <-FALSE
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
  
  
  # Determine dimensions
  Q_X <- ncol(X)
  Q_D <- dim(D)[2] 
  Q_W <- ncol(W)
  Q_A <- dim(A_big)[2] # Can apply to J + 1 portfolios 
  Q_Z <- dim(Z)[2] # Can receive offers from J programs
  N   <- dim(X)[1]
  J   <- dim(Z)[2] # J programs 
  R_draws <- dim(sims_theta)[2] 
  K <- data$K # Number of types for theta
  rm(data)
  gc()
  
  ##############################################################################
  ############################ Extract Parameters ##############################
  ##############################################################################
  
  
  
  
  
  # First parameters correspond to school mean utilities (delta_j)
  l_param <- 1 
  
  if (TRUE){
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
    
    
    # cost parameters (only need fixed cost)
    c_f <- omega[l_param:(l_param + Q_W - 1)]
    l_param <- l_param + Q_W
  } else{
    
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
      -0.3, 0.3,  # ln(sd(theta_1)), ln(sd(theta_2)), #
      0.6932, -1  #alpha_1, mu(theta_1)
    )
    
    delta_j <- omega_choice_0[1:40]
    beta_X <- omega_choice_0[41:52]
    psi_D <- omega_choice_0[53:66]
    c_f <- omega_choice_0[67:79]
  }
  
  # SD of eta_i (goes in costs)
  sigma_c <- exp(omega[l_param])
  l_param <- l_param + 1
  
  # SD of theta_i 
  sigma_theta <- exp(omega[l_param:(l_param + K - 1)])
  l_param <- l_param + K  
  
  
  # alpha_k for population type share
  # We set alpha_K = 0
  if (K>1){
    alpha_k <- c(omega[l_param:(l_param + K - 2)], 0)
    l_param <- l_param + K - 1
  } else if (K==1){
    alpha_k <- 0
  }
  
  # use the soft-max function to convert alphas to p
  exp_alpha_k <- exp(alpha_k)
  sum_exp_alpha_k <- sum(exp_alpha_k)
  p_theta <- exp_alpha_k/sum_exp_alpha_k
  
  # Mu_theta for type k:
  # We normalize mu_K = -(sum^{K-1}_{j=1} p_j mu_j)/p_K
  if (K>1){
    mu_theta <- omega[l_param:(l_param + K - 2)]
    mu_theta <- c(mu_theta, -sum(p_theta[1:(K-1)]*mu_theta)/p_theta[K])
    l_param <- l_param + K - 1
  } else if (K==1){
    mu_theta <- 0
  }
  
  if(posterior){
    log_prior_r <- rep(0, nrow=N, ncol=R_draws)
    theta_r <- matrix(0, nrow = N, ncol = R_draws)
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
    # Extract simulation draws and
    # Transform draws for theta and cost (eta)
    r <- 1
    theta_r <- sweep(outer(as.vector(sims_theta[, r]), sigma_theta), 2, mu_theta, "+")
    eta_ij <- sigma_c * sims_c[,,r]  
    
    ###################### Finish cost calculations ############################
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
    
    
    
    s <- Sys.time()
    
    Emax <- array(0, dim = c(N, (J+1), K))
    
    U_tilde <- array(0, dim=c(N, J, K))
    exp_U_0 <- exp(0)  # Outside option tau is just 0
    exp_U_1 <- exp(U_tilde)  # School options
    
    for (k in 1:K){
      U_tilde[,,k] <- U_bar + matrix(rep(theta_r[,k], J), nrow=N, ncol=J) 
      
      
      # Portfolio 1: Only outside option
      # Only option is outside option
      Emax[, 1, k] <-  0.5772156649  # log(sum) + Euler's constant
      # All other portfolios: 
      for (j in 1:J) {
        # Options are: attend school j or choose outside option + Euler's constant
        Emax[, j+1, k] <- log( exp_U_0 + exp_U_1[, j, k] ) + 0.5772156649
      }
    }
    print(Sys.time()-s)
    
    print(Emax[1:5,1:11,1])
    
    print(Emax[1:5,1:11,2])
    
    s <- Sys.time()
    for (k in 1:K) {
      # Create U_tilde slice
      U_tilde_k <- U_bar + matrix(theta_r[, k], nrow = N, ncol = J)
      U_tilde[,,k] <- U_tilde_k
      
      exp_U_1_k <- exp(U_tilde_k)
      
      # Portfolio 1: Only outside option
      Emax[, 1, k] <- 0.5772156649
      
      # All other portfolios
      for (j in 1:J) {
        Emax[, j + 1, k] <- log1p(exp_U_1_k[, j]) + 0.5772156649
      }
    }
    print(Sys.time()-s)
    
    print(Emax[1:5,1:11,1])
    
    print(Emax[1:5,1:11,2])
    
    
    
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
    if(posterior){
      log_prior_r[,r]  <- -0.5 * (sims_theta_r^2  + sims_c[,r]^2)
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
      dlog_P_A_sd_theta <- matrix(0, nrow = N, ncol = K)
      dlog_P_E_sd_theta <- matrix(0, nrow = N, ncol = K) 
      
      dU <- sigma_theta[k]*sims_theta_r
      dlog_P_A_sd_theta[,k] <- ((1/lambda)*chosenA*dU -  rowSums( (1/lambda) * P_A_short * gA * dU[, rep(1, ncol(gA))], na.rm=TRUE))
      dlog_P_E_sd_theta[,k] <- (ifelse(applied_chose_school, dU, 0) - chosenE*dU)*(1-noEnrollmentStageChoice)
      
      
      
      #################### alpha of theta_k ####################################
      
      # This is one of the terms that'll be later used for the product rule
      dlog_P_A_alpha_K <- matrix(0, nrow = N, ncol = (K-1))
      dlog_P_E_alpha_K <- matrix(0, nrow = N, ncol = (K-1))
      
      if (k == K){
        for (j in 1:(K-1)){
          #print(j)
          dU <- -1* exp_alpha_k[j]*mu_theta[j]
          #print(mu_theta)
          dlog_P_A_alpha_K[,j] <- ((1/lambda)*chosenA*dU)- (rowSums((1/lambda)*gA*dU*P_A_short ,  na.rm=TRUE))
          
          # enrollment stage grad
          dlog_P_E_alpha_K[,j] <- (ifelse(applied_chose_school, 1, 0) - chosenE)*dU*(1-noEnrollmentStageChoice)
        }
      }
      
      #print(dlog_P_A_alpha_K[1:5,])
      
      # print(dlog_P_E_alpha_K[1:5,])
      
      ##################### Mu of theta_i ######################################
      if (K>1){
        dlog_P_A_mu_theta <- matrix(0, nrow = N, ncol = (K-1))
        dlog_P_E_mu_theta <- matrix(0, nrow = N, ncol = (K-1)) 
        
        # Loop over mu_j, take derivative w.r.t mu_j
        for (j in 1:(K-1)){
          
          if (k < K){
            if (j==k){
              # If we're working on type k < K, de/dmu_j is just P_E
              dlog_P_A_mu_theta[,j] <- ((1/lambda)*chosenA)- (rowSums((1/lambda)*gA *P_A_short ,  na.rm=TRUE))
              dlog_P_E_mu_theta[,j] <- (ifelse(applied_chose_school, 1, 0) - chosenE)*(1-noEnrollmentStageChoice)
            }
          } 
          else if (k == K){
            # But if we're working on type k = K, de/dmu_j becomes P_E*(-p_j/p_K)
            de <- (-p_theta[j]/p_theta[K])
            dlog_P_A_mu_theta[,j] <- ((1/lambda)*chosenA* de) - 
              (rowSums((1/lambda)*gA *P_A_short*de, na.rm=TRUE))
            # enrollment stage grad
            dlog_P_E_mu_theta[,j] <- (ifelse(applied_chose_school, 1, 0) - chosenE)*de*(1-noEnrollmentStageChoice)
            
          }
        }
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
      dlog_P_A <- cbind(dlog_P_A_delta, dlog_P_A_X, dlog_P_A_D,dlog_P_A_W,
                        dlog_P_A_sd_eta, dlog_P_A_sd_theta, dlog_P_A_alpha_K, dlog_P_A_mu_theta)
      
      dlog_P_E <- cbind(dlog_P_E_delta, dlog_P_E_X, dlog_P_E_D, dlog_P_E_W, 
                        dlog_P_E_sd_eta, dlog_P_E_sd_theta, dlog_P_E_alpha_K, dlog_P_E_mu_theta)
      
      
      ##################### Compute gradient for draw r #########################
      # Compute total gradient for this simulation draw
      dlikelihood_sim_r <- (dlog_P_A + dlog_P_E) * likelihood_r  
    }
    
    if (grad){
      return(list(likelihood_r=likelihood_r, dlikelihood_sim_r = dlikelihood_sim_r))
    } else {
      return(list(likelihood_r=likelihood_r))
    }
    
  }
  
  # Placeholder for likelihood contribution for i \in I when type = k \in K
  likelihood_i_k <- matrix(0, nrow = N, ncol = K)
  
  # Number of parameters, except for alpha_k
  # Students X N(omega) X Type
  # N(omega) = school mean + X + cost + distance  + sd(theta) + mu(theta) + sd(eta)
  #dlikelihood_i_k <- array(0, dim=c(N, length(omega), K))
  
  dlikelihood_i_k <- array(0, dim=c(N, length(omega), K))
  

  lk_log(sprintf("Computing likelihood when k =  %d...", k))
  res <- mclapply(1:R_draws, function(i){l_and_grad_r(i, k)}, mc.cores=n_cores)
  
  # Create the likelihood matrix (cbind of 1st list elements from each sublist)
  if (grad){
    likelihood_sim <- do.call(cbind, lapply(res, `[[`, 1))
    
    # Create the 3D gradient array from the 2nd list elements
    array_dim <- c(dim(res[[1]][[2]]), length(res))  # N, Q, 30
    dlikelihood_sim <- array(unlist(lapply(res, `[[`, 2)), dim = array_dim)
  } else{
    likelihood_sim <- do.call(cbind, lapply(res, `[[`, 1))
  }
  
  
  # Final progress report
  total_time <- Sys.time() - start_time
  lk_log(sprintf("Completed all %d simulation draws in %.2f minutes", 
                 R_draws, as.numeric(total_time)/60))
  
  # Calculate final likelihood by averaging over simulation draws
  # print(dim(likelihood_sim))
  # print(length(likelihood_sim))
  # print(dim(likelihood_sim))
  likelihood_i <- rowMeans(likelihood_sim, na.rm=TRUE)
  likelihood_i_k[,k] <- likelihood_i
  # print(summary( likelihood_i))
  print(quantile( likelihood_i , probs=c(0.01, 0.05, 0.10, 0.15, 0.20, 0.40, 0.5, 0.60, 0.80, 0.90, 0.95, 0.99), na.rm=TRUE))
  
  # Ensure no zeros (could cause log(0) issues
  threshold_count <- sum( likelihood_i < 1e-6, na.rm=TRUE)
  
  percent_threshold <- round(100*threshold_count/length(likelihood_i), 2)
  lk_log(sprintf("Low likelihood warning: %d observations (%g%%) hit minimum threshold (1e-6)", 
                 threshold_count, percent_threshold))
  lk_log(sprintf("Pre-calculation check: %d NaNs in likelihood_i", sum(is.na(likelihood_i))))
  lk_log(sprintf("Likelihood summary - Min: %.3e, Q1: %.3e, Median: %.3e, Mean: %.3e, Q3: %.3e, Max: %.3e",
                 min(likelihood_i, na.rm=TRUE), quantile(likelihood_i, 0.25, na.rm=TRUE), median(likelihood_i, na.rm=TRUE),
                 mean(likelihood_i, na.rm=TRUE), quantile(likelihood_i, 0.75, na.rm=TRUE), max(likelihood_i, na.rm=TRUE)))
  
  # Calculate final gradient if requested
  if (grad) {
    
    #print(dim(dlikelihood_sim))
    #print(dim(dlikelihood_i_k))
    
    # Average gradient over simulation draws
    dlikelihood_i_k[,,k] <- rowMeans(dlikelihood_sim, dims = 2, na.rm = TRUE)
    #print(summary(dlikelihood_i_k[,,k]))
  }
  
  lk_log("Likelihood function completed")
  
  theta_posterior_mu <- NULL
  theta_posterior_sd <- NULL
  
  # TODO: How do we update posterior means when we have mixture of types?
  if(posterior==TRUE & K == 1){
    # Calculate posterior weights
    log_lik <- log(likelihood_sim)
    log_weights <- log_lik + log_prior
    # Normalize weights
    weights <- exp(log_weights)
    weights <- weights / rowSums(weights)
    
    # Calculate posterior means and variances
    theta_posterior_mu <- rowSums(weights * theta_r)
    theta_posterior_var <- rowSums(weights * (theta_r - theta_posterior_mu)^2)
    theta_posterior_sd <- sqrt(theta_posterior_var)
    
  }
  
  
  # Take the weighted average using p_k as weights
  likelihood_i <- likelihood_i_k %*% matrix(p_theta, nrow=K, ncol = 1)
  
  ## For debugging weighted average by type share
  #print(sum(likelihood_i_k[,1]))
  #print(sum(likelihood_i_k[,2]))
  #print(sum(likelihood_i_k[,3]))
  #print(p_theta)
  #print(sum(likelihood_i))
  #print(sum(-log(likelihood_i)))
  
  
  alpha_k_h <- alpha_k + c(1e-3, 0)
  exp_alpha_h <- exp(alpha_k_h)
  p_theta_h<- exp_alpha_h/sum(exp_alpha_h)
  
  ## For debugging gradients
  #print(alpha_k)
  #print(alpha_k_h)
  #likelihood_i_h <- likelihood_i_k %*% matrix(p_theta_h, nrow=K, ncol = 1)
  #grad_naive<- sum((log(likelihood_i_h) - log(likelihood_i))/1e-3)
  #print(grad_naive)
  #print(sum(log(likelihood_i)))
  #print(sum(log(likelihood_i_h)))
  
  
  if (grad){
    
    
    n_omega <- dim(dlikelihood_i_k)[2]
    dlikelihood_mat <- matrix(dlikelihood_i_k, nrow = (N*n_omega), ncol = K)
    dlikelihood_i <- matrix(dlikelihood_mat %*% as.numeric(p_theta), nrow = N, ncol =n_omega)
    
    # Complete gradients w.r.t, alpha_j 
    # It involves product rule
    # school FE, X, cost, Distance, sd(eta), sd(theta), + 1
    idx <- J+Q_X+Q_W+Q_D+1+K+1
    if (K>1){
      for (j in (1:(K-1))){
        gP <- p_theta[j]*(1-p_theta[j])*likelihood_i_k[,j] 
        
        if (K>2){
          sum_gP_not_j <- rowSums(p_theta[j]*likelihood_i_k[,-j]%*% matrix(p_theta[-j], nrow =(K-1))) 
        } else{
          sum_gP_not_j <- (p_theta[j]*p_theta[-j]*likelihood_i_k[,-j])  #rowSums needed if K>2 only
        }
        
        dlikelihood_i[,idx]<-  gP - sum_gP_not_j + dlikelihood_i_k[,idx,K]*p_theta[K] 
        idx <- idx + 1
      }
      
      
      # Calculate gradient of log likelihood
      likelihood_i_mat <- matrix(likelihood_i, nrow = length(likelihood_i) , ncol = ncol(dlikelihood_i), byrow = FALSE)
      G_i <- dlikelihood_i / likelihood_i_mat # N x No. Parameters
      
    }
    
    # print(quantile(G_i, probs=c(0.05, 0.10, 0.20, 0.40, 0.5, 0.60, 0.80, 0.90, 0.95), na.rm=TRUE))
    # Check for NaN or Inf values
    count_gi_nan <- sum(is.nan(G_i), na.rm = TRUE)
    lk_log(sprintf("Pre-calculation check: %d NaNs in gradient", count_gi_nan  ))
  }
  
  #### Calculate gradients w.r.t alpha_j 
  
  return(list(L = likelihood_i, G = G_i, 
              posterior_mu = theta_posterior_mu, 
              posterior_sd = theta_posterior_sd))
}

