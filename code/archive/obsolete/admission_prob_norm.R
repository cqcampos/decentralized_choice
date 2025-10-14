admission_prob_norm <- function(param,
                                P0,P1,
                                seats,
                                pi_i_init,
                                delta_j, delta_X, psi_D, c_X, 
                                sd_eta, sd_theta_1, sd_theta_2, mu_theta_1, mu_theta_2, p_type_1, p_type_2,
                                preference_model, R=100){
  
  
  ############################ Prepare admission probs #########################
  # param is a J-vector which we will optimize over and correspond to BR admission probs 
  pi <- exp(param) / (1 + exp(param))
  pi <- rep(pi, each=10)
  seats$pi <- pi 
  seats$pi <- seats$pi*seats$share_seats
  # Create endyear x phbao admission probs where columns correspond to schools 
  pi <- seats %>% 
    ungroup() %>% 
    select(endyear, phbao, schoolcode, pi) %>%
    pivot_wider(
      names_from   = schoolcode,
      values_from  = pi,
      names_prefix = "pi_",
      values_fill  = list(pi = 0)   # or NA if you’d rather
    )
  pi_i_init <- as.data.frame(pi_i_init)
  pi_i <- pi_i_init %>%
    left_join(pi, by = c("endyear", "phbao")) %>% arrange(endyear, studentpseudoid)
  
  total_seats <- seats %>%
    group_by(schoolcode) %>%
    summarise(total_seats = sum(total_seats), .groups = 'drop') %>% arrange(schoolcode)
  seats_total <- total_seats$total_seats
  ###################### Calculate application probs #########################
  delta_j_mat <- matrix(delta_j, nrow = N, ncol=length(delta_j), byrow = TRUE)
  delta_X_mat <- X %*% delta_X  
  U_bar <- delta_j_mat + as.numeric(delta_X_mat)

  # Add the distance (D) component (summing over Q_D for each school j)
  D_component <- matrix(0, nrow = N, ncol = J)
  # Takes the sum of each row of D[,,j]* psi_D which is just d_{ij}*X_i*psi_D
  # It is the multiplication by X_i that initially has D[,,,] 
  # After this, U_bar will be almost complete, just missing theta_i 
  for (j in 1:J) {
    D_component[, j] <- rowSums(D[,,j] * matrix(psi_D, nrow = N, ncol = Q_D, byrow = TRUE), na.rm=TRUE)
  }
  U_bar <- U_bar + D_component + theta 

  
  ############################ Application Costs ###############################
  Q_A <- J + 1
  c_X <- matrix(c_X, ncol=1, nrow=length(c_X))
  C_f_bar <- matrix(rep(W %*% c_X, Q_A), ncol = Q_A)
  C_f_bar[,1] <- 0 
  
  # Type-specific application probabilities 
  #q_A_1 <- array(0, dim=c(N, Q_A, R))
  #q_A_2 <- array(0, dim=c(N, Q_A, R))
  
  #p_E_1 <- array(0, dim=c(N, (Q_A-1), R))
  #p_E_2 <- array(0, dim=c(N, (Q_A-1), R))
  pi_i <- pi_i[, 4:ncol(pi_i)]  # Keep only the columns with pi_i for each school
  ncores <- detectCores() - 1
  
  res_list <- mclapply(1:R, function(r) {
#  for(r in 1:R){
    sims_theta <- matrix(rnorm(N * 2, mean = 0, sd = 1), nrow = N, ncol = 2)  # Simulate random coefficients
    sims_c <- array(rnorm(N * Q_A , mean = 0, sd = 1), dim = c(N, Q_A))  # Simulate random costs
    sims_c[,1] <- 0 
    # Add random coefficients to each ubar
    U1 <- U_bar + mu_theta_1 + sd_theta_1 * sims_theta[,1]
    U2 <- U_bar + mu_theta_2 + sd_theta_2 * sims_theta[,2]
    # convert to emax 
    U1e <- pi_i*(log(1 + exp(U1)) ) +   matrix(rep(0.577, times=N*J) , ncol=J, nrow=N) 
    U2e <- pi_i*(log(1 + exp(U2)) +  matrix(rep(0.577, times=N*J) , ncol=J, nrow=N) )
    U1e <- cbind(rep(0.577 ,N), U1e)  # Add outside option utility (0.577 for type 1)
    U2e <- cbind(rep(+ 0.577 ,N), U2e)  # Add outside option utility (0.577 for type 2)
    
    # Add random coefficients to costs 
    C <- exp(C_f_bar + sd_eta * sims_c)
    C[,1] <- 0 
    
    # Calculate the application probabilities for each type
    num <- as.matrix( exp( (U1e - C)/0.1 )  )
    denom <- rowSums2(num, na.rm = TRUE)
    denom <- matrix(rep(denom, times=J+1), ncol = J + 1, nrow=N)
    #q_A_1[,,r] <- num/denom
    q_A_1_r <- num/denom
    num <- as.matrix(exp((U2e - C)/0.1 ))
    denom <- rowSums2(num, na.rm = TRUE)
    denom <- matrix(rep(denom, times=J+1), ncol = J + 1, nrow=N)
    #q_A_2[,,r] <- num/denom
    q_A_2_r <- num/denom
    
    
    # Acceptance probabilities 
    U1 <- U1
    U2 <- U2
    num <- as.matrix(exp(U1))
    denom <- as.matrix((1+ num))
    #p_E_1[,,r] <- num / denom
    p_E_1_r <- num / denom
    num <- as.matrix(exp(U2))
    denom <- as.matrix((1+ num))
    #p_E_2[,,r] <- num / denom
    p_E_2_r <- num / denom
    
    list(
      q_A_1 = q_A_1_r, 
      q_A_2 = q_A_2_r, 
      p_E_1 = p_E_1_r, 
      p_E_2 = p_E_2_r
    )
  }, mc.cores = ncores)
  #}
  q_A_1 <- simplify2array(lapply(res_list, `[[`, "q_A_1"))
  q_A_2 <- simplify2array(lapply(res_list, `[[`, "q_A_2"))
  p_E_1 <- simplify2array(lapply(res_list, `[[`, "p_E_1"))
  p_E_2 <- simplify2array(lapply(res_list, `[[`, "p_E_2"))
  # Average over R simulations
  q_A_1_bar <- rowMeans(q_A_1, dims=2)
  q_A_2_bar <- rowMeans(q_A_2, dims=2)

  p_E_1_bar <- rowMeans(p_E_1, dims=2)
  p_E_2_bar <- rowMeans(p_E_2, dims=2)
  
  # q_A 
  q_A <- p_type_1 * q_A_1_bar + p_type_2 * q_A_2_bar
  q_A <- q_A[, 2:ncol(q_A)]
  q_A_sums <- colSums(q_A)
  print(  sum( q_A_sums )) 
  
  
  # p_E
  p_E <- p_type_1 * p_E_1_bar + p_type_2 * p_E_2_bar
  
  # Now need to sum q_A*p_E*pi_i only Z=1 and A=1
  prod <- q_A * p_E
  E <- E[, 2:ncol(E)]  # Remove the outside option 
  # Calculate the expected number of enrollments
  pred_enroll <- colSums(prod, na.rm = TRUE)
  
  print(pred_enroll)
  init_pi <- exp(param) / (1 + exp(param))
  print(seats_total/pred_enroll)

  pi_br <- pmin(seats_total/pred_enroll, 1) 
  
 
  
  diff <- sqrt(sum((pi_br-init_pi)^2, na.rm=TRUE))
  print(paste0("Current parameter vector: ", paste(param, collapse = ", ")))
  print(paste0("Current predicted BR: ", paste(pi_br, collapse = ", ")))
  print(paste0("Current objective function value: ", diff))
  return(diff)
  
}