convert_to_torch <- function(seats_eq,   
                               stu_type,
                               theta, C_j,
                               delta_j, delta_X, psi_D, c_X, 
                               X, D, W, E, A){
    
    
    # Application end year x phbao indicator matrix for merging pi_ij 
  
    #stu_type_ind_mat[, 1] <- as.integer(stu_type$endyear == 2004 & stu_type$phbao == 0)
    #stu_type_ind_mat[, 2] <- as.integer(stu_type$endyear == 2004 & stu_type$phbao == 1)
    #stu_type_ind_mat[, 3] <- as.integer(stu_type$endyear == 2005 & stu_type$phbao == 0) ...
  
    n <- nrow(stu_type)
    app_years <- sort(unique(stu_type$endyear))
    n_types <- length(app_years) * 2           # number of app years x PHBAO(0,1)
    stu_type_ind_mat <- matrix(0L, n, n_types)
    for (t in app_years){
      for (g in 0:1){
        col_type <-  2*(t-2004) + 1 + g
        #print(paste("Year:", t, "PHBAO:", g, "Type:", col_type ))
        stu_type_ind_mat[, col_type] <-  as.integer(stu_type$endyear == t & stu_type$phbao == g)
      }
    }
    
    
    # seats vector (Assumes that seats are sorted by schoolcode, endyear, phbao in ascending order)
    seats_eq <-  matrix(seats_eq$seats, nrow = n_types)
    
    # Check if CUDA is available
    if (cuda_is_available()) {
      print("Using cuda for objfn computation.")
      device <- torch_device("cuda")
    } else {
      print("Using cpu for objfn computation.")
      device <- torch_device("cpu")
    }
    

    
    # Create a list containing arrays needed for finding best-response pi_ij
    torch_list <- list(
      "seats_eq" = torch_tensor(seats_eq, dtype = torch_float32(), device = device, requires_grad=TRUE),
      "stu_type" = torch_tensor(stu_type_ind_mat, dtype = torch_float32(), device = device, requires_grad=TRUE),
      "theta"    = torch_tensor(theta, dtype = torch_float64(), device = device, requires_grad=TRUE),
      "C_j"      = torch_tensor(C_j, dtype = torch_float64(), device = device, requires_grad=TRUE),
      "X"        = torch_tensor(X, dtype = torch_float64(), device = device, requires_grad=TRUE),
      "D"        = torch_tensor(D, dtype = torch_float64(), device = device, requires_grad=TRUE),
      "W"        = torch_tensor(W, dtype = torch_float64(), device = device, requires_grad=TRUE),
      "E"        = torch_tensor(E, dtype = torch_float64(), device = device, requires_grad=TRUE),
      "A"        = torch_tensor(A, dtype = torch_float64(), device = device, requires_grad=TRUE),
      "delta_j"  = torch_tensor(delta_j, dtype = torch_float64(), device = device, requires_grad=TRUE),
      "delta_X"  = torch_tensor(delta_X, dtype = torch_float64(), device = device, requires_grad=TRUE),
      "psi_D"    = torch_tensor(psi_D, dtype = torch_float64(), device = device, requires_grad=TRUE),
      "c_X"      = torch_tensor(c_X, dtype = torch_float64(), device = device, requires_grad=TRUE)
    )
    
    
    
    
    
  }
  
# Takes param_tens (tensor)  as an input
pi_norm_torch_grad <- function(param_tens, inp, zero_seat_bool, preference_model, eta){
  
    # s <- Sys.time()
    ############# Merge endyear x phbao specific pi_ij to students ###############
    # param is a (JXNTYPES)-vector which we will optimize over and correspond to BR admission probs 
    # 10 = number of endyear (5) x phbao (2)
    
   
    # param[1] : 1st school, year 2004, non-phbao
    # param[2] : 1st school, year 2004, phbao
    # param[3] : 1st school, year 2005, non-phbao
    # Create a student-level data, where jth column corresponds to pi_g(i)j 
    
    # K = number of end year x phbao 
    K <- inp$stu_type$shape[2]
    J <- length(zero_seat_bool)/K
  
    
    # student-by-school matrix of current pi
    # Rows correspond to student, columns correspond school's admission pr, based on i's endyear x phbao
    # param_tens <- torch_exp(param_tens)/(1+exp(param_tens))
      
    # Initialize array of Types x J, set grad = False
    # 
    pi_types_tens <- torch_zeros(c(length(zero_seat_bool)), 
                                 requires_grad=TRUE, 
                                 device = inp$X$device)
    
    #print(length(zero_seat_bool))
    #print( list(torch_tensor(which(!zero_seat_bool),dtype=torch_long())))
    pi_types_tens <- pi_types_tens$index_put(
                          list(torch_tensor(which(!zero_seat_bool),dtype=torch_long())), 
                          param_tens)
 
    
    pi_types_tens <- pi_types_tens$reshape(c(J, K))$t()
    

    pi_i <- torch_matmul(inp$stu_type, pi_types_tens)   # (N x TYPES) %*% (TYPES x J) -> (N x J)
    
    
  
    ################## Calculate P_A (Prob. choosing each portfolio) #####################
    
    # Calculate utility of each portfolio during the app stage
    # u_bar is N by J matrix, where u_bar[i,j] is i's app-stage utility of choosing school j
    
    # CC Change
    #u_bar <- inp$delta_j$unsqueeze(1) +  # school-specific mean utility (None by J)
    #  torch_matmul(inp$X, inp$delta_X)$unsqueeze(2) + # decision maker mean utility (N by None)
    #  torch_einsum("nkj,k->nj", list(inp$D, inp$psi_D)) + # distance effects (N by J)
    #  inp$theta$unsqueeze(2)                              # Individual r.c. (N by none)
    u_bar <- inp$delta_j$unsqueeze(1) +  # school-specific mean utility (None by J)
      torch_matmul(inp$X, inp$delta_X)$unsqueeze(2)  # decision maker mean utility (N by None)
    # Add distance effects if not eliminating travel costs 
    if(preference_model !=5 & preference_model!=6 & preference_model != 9){
      u_bar <- u_bar + torch_einsum("nkj,k->nj", list(inp$D, inp$psi_D)) # distance effects (N by J)
    }
    u_bar <- u_bar + inp$theta$unsqueeze(2) # Individual r.c. (N by none) 
    # Calculate net expected utility max
    # Define helper function for calculating cost
    cost_app <- function(preference_model, inp, eta_spec){
      if (preference_model == 1){
        C <- 0
      } else{
        
        # length n vector of mean cost
        c_f_bar <- torch_einsum("nk,k -> n", list(inp$W, inp$c_X))
        
        if (eta_spec == "eta_in"){
          return(torch_exp(c_f_bar + inp$C_j)$unsqueeze(2) ) # cost shock
        } else if (eta_spec == "eta_out"){ # Ryan: Added!
          return(torch_exp(c_f_bar)$unsqueeze(2) + inp$C_j)
        } else{
          print("Error from cost_app(): Please set 'eta_spec' to 'eta_in' or 'eta_out'!!!")
        }
      }
    }
    
    
    
    # note: 1-pi_ij not shown in the computation, as it's cancelled out 
    lambda <- 0.05
    euler_cons <- .577
    
    emax <- pi_i*torch_log(1+torch_exp(u_bar)) + euler_cons
    
    C    <-  cost_app(preference_model,
                      inp,
                      eta)
    

    net_util <- emax - C
      
    # Add a column for outside option (not applying)
    util_out <- torch_full(c(net_util$size()[1], 1), # N x 1
                           euler_cons)$to(device = inp$X$device)        # u(outside)
    
    net_util <- torch_cat(list(util_out, net_util), dim = 2)
  
    num <- torch_exp(net_util/lambda)
    P_A <- num/(num$sum(dim=2, keepdim=TRUE))
    
    
    ################################# P_E ########################################

    num <-  torch_exp(u_bar)
    P_E <- num/(1+num)

    
    # Create endyear x phbao expected enrollment counts where columns correspond to schools 
    # stu_type <- N by K matrix
    # P_A[,2:41]*P_E* pi_i <- N by J matrix
    # pred_enroll[k, j] <- sum_n\in k p_ae_pi[n, j]
    # n is in k if nstu_type[n, k] == 1
    
  # p_ae_pi <- (P_A[,2:(P_A$size()[2])]* P_E *pi_i)
    p_a_e <-  (P_A[,2:(P_A$size()[2])]* P_E) # predicted enrollment | offer
    stu_type_f <- inp$stu_type$to(dtype = p_a_e$dtype)
    pred_enroll <- stu_type_f$transpose(1, 2)$matmul(p_a_e)  # (K, J)
    pred_app <- stu_type_f$transpose(1,2)$matmul(P_A[,2:(P_A$size()[2])]) # (K, J)
 
    
    # Calculate pi Best Response
    eps <- torch_tensor(1e-12, dtype = pred_enroll$dtype, device = pred_enroll$device, requires_grad=TRUE)
    
    # Safe denominator
    denom <- torch_where(pred_enroll$gt(0), pred_enroll, eps)
    
    # Compute ratio only where denom > 0; 0 otherwise (avoids NaN/Inf in the graph)
    pi_raw <- torch_where(
      pred_enroll$gt(0),
      inp$seats_eq / denom,
      torch_zeros_like(pred_enroll)
    )
    
    # Apply your rules WITHOUT in-place indexing
    pi_br <- torch_clamp(pi_raw, min = 0, max = 1)
    
    # Force to 1 where seats>0 & pred_enroll==0
    mask_pos_zero <- pred_enroll$eq(0) & inp$seats_eq$gt(0)
    pi_br <- torch_where(mask_pos_zero, torch_ones_like(pi_br), pi_br)
    
    # Force to 0 where seats==0
    mask_zero_seat <- inp$seats_eq$eq(0)
    pi_br <- torch_where(mask_zero_seat, torch_zeros_like(pi_br), pi_br)
    
  
    
    #pi_br[pi_br > 1] <- 1
    #pi_br[pred_enroll == 0 & inp$seats_eq > 0] <- 1
    #pi_br[inp$seats_eq ==0] <- 0
    
    

    
    # Objfn value:
    # diff <- (pi_br - param_tens)
  
    diff <- (pi_br$transpose(1,2)$reshape(length(zero_seat_bool))[which(!zero_seat_bool)] - param_tens)

    #print("=====================================================")
    
    #print("Predicted Enrollments:")
    #print(pred_enroll, n=-1)
    
    #print("Current pi:")
    #print(param_tens, n = -1)
    
    #print("Best response pi:")
    #print(pi_br, n=-1)
  
    
    # diff <- (diff$norm(p = 2))^2
    # cat("Current objective function value:\n", diff, "\n")
    
    
    
    # print(paste0("Time elapsed:", Sys.time()-s))
    #return(2*as.numeric(diff))
    # diff <- (diff$pow(2)$sum())
    
    
    diff <- (diff$pow(2)$sum())
    return(diff)
  
    
  }
  

pi_norm_torch <- function(param, inp, zero_seat_bool, preference_model, eta, iter){

  
  # s <- Sys.time()
  ############# Merge endyear x phbao specific pi_ij to students ###############
  
  # param is a (JX10)-vector which we will optimize over and correspond to BR admission probs 
  # 10 = number of endyear (5) x phbao (2)
  
  # param[1] : 1st school, year 2004, non-phbao
  # param[2] : 1st school, year 2004, phbao
  # param[3] : 1st school, year 2005, non-phbao
  # Create a student-level data, where jth column corresponds to pi_g(i)j 
  
  # Note: The gradients are undefined at 0 for schools with 0 seats
  #       For those schools,don't estimate pi_br. Plug in pr = 0  for those schools
  param_full <- numeric(length(zero_seat_bool))
  param_full[zero_seat_bool] <- 0 
  param_full[!zero_seat_bool] <- param
  param <- param_full
  

  
  # K = number of end year x phbao 
  K <- ncol(inp$stu_type)
  stopifnot(length(param) %% K == 0L)
  J <- length(param) / K
  
  # reshape param so column j holds the 10 pi_i for school j
  #param <- exp(param)/(1+exp(param))
  param_mat <- matrix(param, nrow = K)   # TYPES x J, column-wise fill
  param_tens <- torch_tensor(param_mat, dtype = torch_float32())$to(device = inp$X$device)
  
  
  # student-by-school matrix of current pi
  # Rows correspond to student, columns correspond school's admission pr, based on i's endyear x phbao
  # Construct a worker-specific logfile name using PID
  pi_i <- torch_matmul(inp$stu_type, param_tens)   # (N x TYPES) %*% (TYPES x J) -> (N x J)
  
  
  ################## Calculate P_A (Prob. choosing each portfolio) #####################
  
  # Calculate utility of each portfolio during the app stage
  # u_bar is N by J matrix, where u_bar[i,j] is i's app-stage utility of choosing school j
  
  #u_bar <- inp$delta_j$unsqueeze(1) +  # school-specific mean utility
  #  torch_matmul(inp$X, inp$delta_X)$unsqueeze(2) + # decision maker mean utility
  #  torch_einsum("nkj,k->nj", list(inp$D, inp$psi_D)) + # distance effects
  #  inp$theta$unsqueeze(2)
  
  u_bar <- inp$delta_j$unsqueeze(1) +  
    torch_matmul(inp$X, inp$delta_X)$unsqueeze(2) + 
    inp$theta$unsqueeze(2)

  if (preference_model != 5 & preference_model != 6 & preference_model != 9) {
    distance_component <- torch_einsum("nkj,k->nj", list(inp$D, inp$psi_D))
    u_bar <- u_bar + distance_component
  }
  


  # Calculate net expected utility max
  # Define helper function for calculating cost
  cost_app <- function(preference_model, inp, eta_spec){
    if (preference_model == 1){
      C <- 0
    } else{
      
      # length n vector of mean cost
      c_f_bar <- torch_einsum("nk,k -> n", list(inp$W, inp$c_X))
      
      if (eta_spec == "eta_in"){
        return(torch_exp(c_f_bar)$unsqueeze(2) + inp$C_j) # cost shock
      } else if (eta_spec == "eta_out"){ # Ryan: Added!
        return(torch_exp(c_f_bar)$unsqueeze(2) + inp$C_j)
      } else{
        print("Error from cost_app(): Please set 'eta_spec' to 'eta_in' or 'eta_out'!!!")
      }
    }
  }
  
  
  
  # note: 1-pi_ij not shown in the computation, as it's cancelled out 
  lambda <- 0.05
  euler_cons <- .577
  emax <- pi_i*torch_log(1+torch_exp(u_bar)) + euler_cons
  
  C    <-  cost_app(preference_model,
                    inp,
                    eta)
  
  
  net_util <- emax - C
  
  # Add a column for outside option (not applying)
  util_out <- torch_full(c(net_util$size()[1], 1), # N x 1
                         euler_cons)$to(device = inp$X$device)        # u(outside)
  
  net_util <- torch_cat(list(util_out, net_util), dim = 2)
  
  num <- torch_exp(net_util/lambda)
  P_A <- num/(num$sum(dim=2, keepdim=TRUE))
  
  
  ################################# P_E ########################################
  
  num <-  torch_exp(u_bar)
  P_E <- num/(1+num)
  
  
  # Create endyear x phbao expected enrollment counts where columns correspond to schools 
  # stu_type <- N by K matrix
  # P_A[,2:41]*P_E* pi_i <- N by J matrix
  # pred_enroll[k, j] <- sum_n\in k p_ae_pi[n, j]
  # n is in k if nstu_type[n, k] == 1
  
  #predicted enrollment conditional on offer
  p_a_e <- (P_A[,2:(P_A$size()[2])]* P_E) 
  stu_type_f <- inp$stu_type$to(dtype = p_a_e$dtype)
  
  pred_enroll <- stu_type_f$transpose(1, 2)$matmul(p_a_e)  # (K, J)
  

  # Calculate pi Best Response
  diff <- inp$seats_eq[inp$seats_eq>0]
  pi_br <- inp$seats_eq/pred_enroll
  pi_br[pi_br > 1] <- 1
  pi_br[pred_enroll == 0 & inp$seats_eq > 0] <- 1
  pi_br[inp$seats_eq ==0] <- 0
  

  # Objfn value:
  pi_br_vec <- pi_br$transpose(1,2)$reshape(length(zero_seat_bool))[which(!zero_seat_bool)] 
  curr_pi_vec <- param_tens$transpose(1,2)$reshape(length(zero_seat_bool))[which(!zero_seat_bool)]
  diff <- (pi_br_vec - curr_pi_vec)
  
  
  
  #print("=====================================================")
  
#  pred_enroll_vec <- pred_enroll$transpose(1,2)$reshape(length(zero_seat_bool))[which(!zero_seat_bool)]
  #seats_mat       <- seats_eq[which(!zero_seat_bool),]
#  debug_df <- data.frame("pi_br" =  as.numeric(pi_br_vec), 
#                         "curr_pi" = as.numeric(curr_pi_vec), 
#                         "pred_enroll" = as.numeric(pred_enroll_vec))
#  debug_df <- cbind(seats_mat, debug_df) 
#  print(debug_df, n = 50)
  
  
  diff <- (diff$pow(2)$sum())$item() 
  cat("Round", iter, "Current objective function value:", diff, "\n")
  
  # print(paste0("Time elapsed:", Sys.time()-s))
  # print(diff)

  return(diff)
  
  
}

