################################################################################
# Script: calculate_eq.R 
# Author: Chris Campos 
#
# Estimates equilibrium admission probs 
# 
# Dependencies:

################################################################################
calculate_eq <- function(param, eta_spec, seats_eq, pi_i, delta_j, delta_X, psi_D, c_X,
                         theta, C_j, 
                         X, D, W, A, E, 
                         stu_type, round,
                         preference_model, use_torch_for_eq = FALSE
                         ){
  
  options_fmin <- list(maxit = 1000,  factr=4.5e10)
  

  if (use_torch_for_eq){
    tensor_list <- convert_to_torch(seats_eq,   
                                    stu_type,
                                    theta, C_j,
                                    delta_j, delta_X, psi_D, c_X, 
                                    X, D, W, E, A)
    
    # Check if CUDA is available
    if (cuda_is_available()) {
      print("Using cuda for objfn computation.")
      device <- torch_device("cuda")
    } else {
      print("Using cpu for objfn computation.")
      device <- torch_device("cpu")
    }
    
    
    # Helper function for calculating gradients using autodiff
    fn_grad <- function(param, tensor_list, zero_seat_bool, preference_model, eta, device, iter){
      
      param_tens <- torch_tensor(param, dtype = torch_float32(), requires_grad=TRUE, device = device)
      
  
      
      objfn_tens <- pi_norm_torch_grad(param_tens,tensor_list, zero_seat_bool, preference_model,eta)
    
      cat("Round", iter, "Objfn returned from fn_grad", objfn_tens$item(), "\n")
      
      objfn_tens$backward()
      
      grad_vec  <- param_tens$grad        # flatten to 400 (1D tensor)
      grad_vecR <- as.numeric(grad_vec)       # convert to base R numeric vector
      param_tens$detach()

      grad_norm <- sqrt(sum(grad_vecR^2, na.rm = TRUE))
      
    # Warn if NaN values present
        if (any(is.nan(grad_vecR))) {
        cat("NaN detected in gradient vector!")
          print(grad_vecR)
      }
      
      cat("Round", iter, "Current Gradient Norm:", grad_norm, "\n")
      
      return (grad_vecR)
    }
    
    
    # Note: The gradients are undefined at 0 for schools with 0 seats
    #       For those schools,don't estimate pi_br 
    # Only solve for schools with seats
    
    zero_seat_bool <- seats_eq$seats==0
    #param <- runif(n=length(param), min = 0.1, max =1)
    #param <- runif(n=length(zero_seat_bool))
    param_solve <- param[which(!zero_seat_bool)]
    #print("First 10 values of pi we're solving")
    #print(param_solve[1:10])
    pi_br_optim <- optim(param_solve,
                         fn = function(x) { pi_norm_torch(x, 
                                                          tensor_list,
                                                          zero_seat_bool,
                                                          preference_model,
                                                          eta,round)
                         },
                         gr = function(x){
                           fn_grad(x, 
                                   tensor_list,
                                   zero_seat_bool,
                                   preference_model,
                                   eta, device, round)
                         },
                         method = "L-BFGS-B",
                         lower = 0,
                         upper = 1,
                         control = options_fmin)
    
    # Save so that next iteration starts from a better point
    pi_br_solved <- numeric(length(zero_seat_bool))
    pi_br_solved[zero_seat_bool] <- 0 
    pi_br_solved[!zero_seat_bool] <- pi_br_optim$par
    
    write.csv(pi_br_solved, paste0(dir, "/estimates/", eta_spec, "_K", K, "_init_param_preference_model_",preference_model,"_", last_yr, ".csv"), row.names = FALSE)
    
    

  } else{
  
    pi_br_optim <- optim(param,
                         fn = function(x) { pi_norm(x, 
                                                    eta_spec = eta_spec,
                                                    seats_eq=seats_eq,
                                                    stu_type = stu_type,
                                                    delta_j = delta_j, delta_X = delta_X,
                                                    psi_D = psi_D, c_X = c_X,
                                                    theta = theta, C_j = C_j,
                                                    X = X, D = D, W = W, E = E, A = A,
                                                    preference_model = preference_model)
                         },
                         method = "L-BFGS-B",
                         control = options_fmin)
    
    pi_br_solved <- pi_br_optim$par
  }

  
  return(pi_br_solved)
  
  
}
