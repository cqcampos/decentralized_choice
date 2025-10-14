pi_norm <- function(param,
                    eta_spec,
                    seats_eq,
                    stu_type,
                    delta_j, delta_X, psi_D, c_X, 
                    theta, C_j, 
                    X, D, W, E, A, 
                    preference_model){
  
  
  ############################ Prepare admission probs #########################
  # param is a (JX10)-vector which we will optimize over and correspond to BR admission probs 
  # 10 = number of endyear (5) x phbao (2)
  seats_eq$pi <- param
  
  # Create endyear x phbao admission probs where columns correspond to schools 
  pi_mat <- seats_eq %>% 
    ungroup() %>% 
    select(endyear, phbao, schoolcode, pi) %>%
    pivot_wider(
      names_from   = schoolcode,
      values_from  = pi,
      names_prefix = "pi_",
      values_fill  = list(pi = 0)   # or NA if you’d rather
    )

  # Merge the school x endyear x phbao probability to students
  pi_i <- stu_type %>%
    left_join(pi_mat, by = c("endyear", "phbao")) 
  

  ###################### Calculate application probs #########################
  delta_j_mat <- matrix(delta_j, nrow = N, ncol=length(delta_j), byrow = TRUE)
  delta_X_mat <- X %*% delta_X  
  U_bar <- delta_j_mat + as.numeric(delta_X_mat)

  # Add the distance (D) component (summing over Q_D for each school j)
  #D_component <- matrix(0, nrow = N, ncol = J)
  # Takes the sum of each row of D[,,j]* psi_D which is just d_{ij}*X_i*psi_D
  # It is the multiplication by X_i that initially has D[,,,] 
  # After this, U_bar will be almost complete, just missing theta_i 
  #for (j in 1:J) {
  #  D_component[, j] <- rowSums(D[,,j] * matrix(psi_D, nrow = N, ncol = Q_D, byrow = TRUE), na.rm=TRUE)
  #}
  
  # D: N x Q_D x J, psi_D: length Q_D
  # Verified that does same as above
  Dj <- aperm(D, c(1, 3, 2))        # N x J x Q_D  (put Q last)
  dim(Dj) <- c(N * J, Q_D)          # (N*J) x Q_D
  
  # replicate na.rm=TRUE from rowSums:
  Dj[is.na(Dj)] <- 0
  D_component <- matrix(Dj %*% psi_D, nrow = N, ncol = J)
  
  U_bar <- U_bar + D_component + theta 
  
  
  ############################ Application Costs ###############################
  Q_A <- J + 1
  if(preference_model!=1){
    c_X <- matrix(c_X, ncol=1, nrow=length(c_X))
    C_f_bar <- matrix(rep(W %*% c_X, Q_A), ncol = Q_A)
    C_f_bar[,1] <- 0 
    
    C_j <- cbind(rep(0, times=N), C_j) # Add the outside option cost of 0
    #print(dim(C_j))
    #print(dim(C_f_bar))
    # Note from Chris: Should add a condition for eta_in vs eta_out
    if (eta_spec == "eta_in"){
      C <- exp(C_f_bar + C_j) # cost shock
    } else if (eta_spec == "eta_out"){ # Ryan: Added!
      C <- exp(C_f_bar) + C_j
    } else{
      print("PLEASE SET eta_spec TO eta_in OR eta_out!!!!")
    }
    C[,1] <- 0
  }else{
    C <- 0
  }
  
  ################################# P_A ########################################
  # Calculate net expected utility 
  pi_i <- pi_i[, 4:ncol(pi_i)]
  lambda <- 0.05
  Ue <-  (pi_i*(log(1 + exp(U_bar)) ) +   matrix(rep(0.577, times=N*J) , ncol=J, nrow=N))
  Ue <- cbind(rep(0.577 ,N), Ue)  # Add outside option utility (0.577)
  num <- as.matrix( exp( (Ue - C)/lambda )  )
  denom <- rowSums2(num, na.rm = TRUE)
  denom <- matrix(rep(denom, times=J+1), ncol = J + 1, nrow=N)
  P_A <- num/denom
  
  ################################# P_E ########################################
  num <- as.matrix(exp(U_bar))
  denom <- as.matrix((1+ num))
  P_E <- num / denom

  # Now need to sum q_A*p_E*pi_i only Z=1 and A=1
  P_A <- P_A[, 2:ncol(P_A)]
  
  # Create endyear x phbao expected enrollment counts where 4th ~ last columns correspond to schools 
  stu_type <- cbind(stu_type, (P_A * P_E *pi_i) )
  
  pred_enroll <- stu_type %>% 
    group_by(endyear, phbao) %>%
    summarise(across(starts_with("pi_"), sum, na.rm = TRUE), 
              .groups = "drop")
  

  # Flatten expected enrollments into 1D array
  pred_enroll <- as.vector(as.matrix(pred_enroll[,3:ncol(pred_enroll)]))

  
  # Calculate pi Best Response
  seats_total <- seats_eq$seats
  pi_br <- pmin(as.numeric(seats_total/pred_enroll), 1) 
  
  
  # Objfn value:
  diff <- sqrt(sum((pi_br-param)^2, na.rm=TRUE))
  
  
  # Print out numbers calculated during the optimization for the first 20 school x year x phbao
  seats_eq$pi_curr <- param 
  seats_eq$pi_br   <- pi_br 
  seats_eq$pred_enr <- pred_enroll
  print(seats_eq[1:20, c("endyear", "phbao", "schoolcode", "seats", "pred_enr", "pi_br", "pi_curr")])
  
  cat("Current objective function value:\n", diff, "\n")
  
  
  return(diff)
  
}