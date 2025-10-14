# Function for evaluating model's accuracy
sim_a_e <- function(endyear, low_score_school, mean_school_scores, 
                                phbao, asian, white, 
                                yearfe, blockfe, outcome_fes, 
                                alpha_j_math, alpha_0_X_math, alpha_m_X_math, alpha_0_math, alpha_m_math,
                                alpha_j_ela, alpha_0_X_ela, alpha_m_X_ela, alpha_0_ela, alpha_m_ela,
                                delta_j, delta_X, psi_D, c_X, 
                                sd_eta, sd_theta_1, sd_theta_2, sd_theta_3,
                                mu_theta_1, mu_theta_2, mu_theta_3,
                                p_type_1, p_type_2, p_type_3,
                                params, exog,eta ) {
  
  # 1. Unpack exogenous data
  X <- exog$X            # N×Q_X
  D <- exog$D            # N×Q_D×J
  W <- exog$W            # N×Q_W
  seats <- exog$seats 
  
  
  N <- nrow(D)
  J <- dim(D)[3]
  Q_X <- ncol(X) 
  Q_D <- dim(D)[2] 
  Q_W <- ncol(W)
  
  # 2. Unpack some details 
  Amat            <- params$Amat  # NA×J
  pi_i            <- params$pi_i # equilibrium π_i
  K               <- params$K  # number of types
  homog           <- params$homog  # homogeneous flag
  preference_model<- params$preference_model
  Q_A <- ncol(Amat) 
  Q_Z <- J + 1 
  
  
  # 3. Simulate shocks 
  # post-lottery shock
  xi   <- -log(-log(matrix(runif(N*(J+1)), N, J+1)))
  # cost shock
  C_j  <- matrix(rnorm(N*J, 0, sd_eta), N, J)
  
  # 4. Simulate latent types and theta_i
  p_theta <- c(p_type_1, p_type_2, p_type_3)
  type_draw <- sample(1:3, size = N, prob = p_theta, replace=TRUE) 
  if(K==1){ 
    theta_choice <- mu_theta_1 + sd_theta_1 * rnorm(N) 
    theta <- theta_choice 
  }else if(K==2){
    theta_1 <- mu_theta_1 + sd_theta_1 * rnorm(N) 
    theta_2 <- mu_theta_2 + sd_theta_2 * rnorm(N) 
    theta_choice <- ifelse(type_draw == 1, theta_1, theta_2)
    theta <- theta_choice 
    rm(theta_1, theta_2) 
  }else if(K==3){
    theta_1 <- mu_theta_1 + sd_theta_1 * rnorm(N)
    theta_2 <- mu_theta_2 + sd_theta_2 * rnorm(N)
    theta_3 <- mu_theta_3 + sd_theta_3 * rnorm(N)
    theta_choice <- ifelse(type_draw == 1, theta_1, ifelse(type_draw == 2, theta_2, theta_3)) 
    theta <- theta_choice 
    rm(theta_1, theta_2, theta_3)
  }
  
  
  # 5. Compute utilities of attending each school
  # Construct mean utilities (U_bar) 
  # We have a constant psi_j (group_matrix) so let's make a matrix of that constant 
  delta_j_mat <- matrix(delta_j, nrow = N, ncol=length(delta_j), byrow = TRUE)
  delta_X_mat <- X %*% delta_X  
  
  # let's add them together and that's our U_bar 
  U_bar <- sweep(delta_j_mat, 1, delta_X_mat, FUN = "+")
  
  # Add the distance (D) component (summing over Q_D for each school j)
  D_component <- matrix(0, nrow = N, ncol = J)
  # Takes the sum of each row of D[,,j]* psi_D which is just d_{ij}*X_i*psi_D
  # It is the multiplication by X_i that initially has D[,,,] have the extra 
  # dimension 
  for (j in 1:J) {
    D_component[, j] <- rowSums(D[,,j] * matrix(psi_D, nrow = N, ncol = Q_D, byrow = TRUE))
  }
  U_bar <- U_bar + D_component
  

  U_tilde <- U_bar  + matrix(theta_choice, N, J) 
  
  #6. Emax calculations
  exp_U0      <- exp(0)
  exp_U_1  <- exp(U_tilde)
  Emax <- matrix(0, nrow = N, ncol = J+1)
  
  # Portfolio 1: No offers received
  # Only option is outside option
  Emax[, 1] <- log(exp_U0) + 0.577  # log(sum) + Euler's constant
  
  # For each school j (portfolio j+1 = offer from school j only)
  for (j in 1:J) {
    # Options are: attend school j or choose outside option
    Emax[, j+1] <- log(exp_U0 + exp_U_1[, j]) + 0.577
  }
  
  # 7. Compute application costs C
  c_X <- matrix(c_X, ncol=1, nrow=length(c_X))
  C_f_bar <- matrix( rep(exp(W %*% c_X), Q_A), ncol = Q_A) 
  if(eta =="eta_in"){
    C <- C_f_bar * exp(C_j) 
  }else{
    C <- C_f_bar + C_j 
  }
  
  # Set application costs to zero if preference_model==1 
  if(preference_model==1){ 
    C <-  matrix(0, nrow = N, ncol = Q_A)
  }
  
  # 8. Compute expected utility of each portfolio
  V <- matrix(0, nrow = N, ncol = Q_A+1)
  # Portfolio 1: Not applying anywhere
  # Expected utility is just the utility of the outside option
  V[, 1] <- Emax[, 1]  # No applications → only get the "no offers" outcome

  for (j in 1:J) {
    V[, j+1] <- pi_i[, j] * Emax[, j+1] + (1 - pi_i[, j]) * Emax[, 1] - C[, j]
  }
  
  # 9. Application choice: compute probability of each portfolio being chosen
  lambda <- .05
  num <- as.matrix(exp(V/lambda))
  denom <- rowSums2(num, na.rm = TRUE)
  denom <- matrix(rep(denom, times=J+1), ncol = J+1, nrow=N)
  P_A <- num/denom
  
  # 10. Enrollment choice: compute probability of enrolling into a school | offer
  num <- as.matrix(exp(U_tilde))
  denom <- as.matrix(1+num)
  P_E <- num / denom
  
  
  # 11. Calculate stats for model accuracy evaluations
  n_apps <- colSums(P_A)  # E[N application for each school]
  n_offers <- colSums(P_A[,(2:(J+1))] * pi_i) # E[N offers given to the applicants]
  n_enroll <- colSums(P_A[,(2:(J+1))] * pi_i * P_E) # E[N enrollments]
  
  n_apps_phbao <- colSums(P_A[which(phbao==1),])  # E[N application for each school]
  n_apps_asian <- colSums(P_A[which(asian==1),])  # E[N application for each school]
  n_apps_white <- colSums(P_A[which(white==1),])  # E[N application for each school]
  n_apps_low_score_school <- colSums(P_A[which(low_score_school==1),])  # E[N application for each school]
  n_apps_high_score_school <- colSums(P_A[which(low_score_school==0),])  # E[N application for each school]
  
  n_offers_phbao <- colSums(P_A[which(phbao==1),(2:(J+1))] * pi_i[which(phbao==1),]) # E[N offers given to the applicants]
  n_offers_asian <- colSums(P_A[which(asian==1),(2:(J+1))] * pi_i[which(asian==1),]) # E[N offers given to the applicants]
  n_offers_white <- colSums(P_A[which(white==1),(2:(J+1))] * pi_i[which(white==1),]) # E[N offers given to the applicants]
  n_offers_low_score_school <- colSums(P_A[which(low_score_school==1),(2:(J+1))] * pi_i[which(low_score_school==1),]) # E[N offers given to the applicants]
  n_offers_high_score_school <- colSums(P_A[which(low_score_school==0),(2:(J+1))] * pi_i[which(low_score_school==0),]) # E[N offers given to the applicants]
  
  n_enroll_phbao <- colSums(P_A[which(phbao==1),(2:(J+1))] * pi_i[which(phbao==1),] * P_E[which(phbao==1),]) # E[N enrollments]
  n_enroll_asian <- colSums(P_A[which(asian==1),(2:(J+1))] * pi_i[which(asian==1),] * P_E[which(asian==1),]) # E[N enrollments]
  n_enroll_white <- colSums(P_A[which(white==1),(2:(J+1))] * pi_i[which(white==1),] * P_E[which(white==1),]) # E[N enrollments]
  n_enroll_low_score_school <- colSums(P_A[which(low_score_school==1),(2:(J+1))] * pi_i[which(low_score_school==1),] * P_E[which(low_score_school==1),]) # E[N enrollments]
  n_enroll_high_score_school <- colSums(P_A[which(low_score_school==0),(2:(J+1))] * pi_i[which(low_score_school==0),] * P_E[which(low_score_school==0),]) # E[N enrollments]
  
  return(list(n_apps = n_apps, n_offers = n_offers, n_enroll = n_enroll, 
              P_E, U_bar = U_bar, U_bar_md = (U_bar-delta_j_mat),
              n_apps_phbao = n_apps_phbao, n_apps_asian = n_apps_asian, n_apps_white = n_apps_white,
              n_apps_low_score_school = n_apps_low_score_school, n_apps_high_score_school = n_apps_high_score_school,
              n_offers_phbao = n_offers_phbao, n_offers_asian = n_offers_asian, n_offers_white = n_offers_white,
              n_offers_low_score_school = n_offers_low_score_school, n_offers_high_score_school = n_offers_high_score_school,
              n_enroll_phbao = n_enroll_phbao, n_enroll_asian = n_enroll_asian, n_enroll_white = n_enroll_white,
              n_enroll_low_score_school = n_enroll_low_score_school, n_enroll_high_score_school = n_enroll_high_score_school
              )
         
         )
}
