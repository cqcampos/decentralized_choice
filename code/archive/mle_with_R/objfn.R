###################### Helper Function - objfn ###############################
objfn <- function(omega, data, grad = FALSE) {
  # Start timer
  start_time <- Sys.time()
  
  # Unpack data
  Nblocks <- data$Nblocks
  A       <- data$A
  A_dum   <- data$A_dum
  Z       <- data$Z
  E       <- data$E
  X       <- data$X
  D       <- data$D
  W       <- data$W
  A_big   <- data$A_big
  Z_big   <- data$Z_big
  pi_i    <- data$pi_i
  lambda  <- data$lambda
  homog   <- data$homog
  sims_theta <- data$sims_theta
  sims_c     <- data$sims_c
  cl        <- data$cl # Extract cluster object if it exists
  log_file  <- data$log_file  # Extract log file path
  

  # Log function
  obj_log <- function(msg) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    formatted_msg <- paste0("[objfn - ", timestamp, "] ", msg)
    # We use cat with file argument to avoid locking issues
    cat(formatted_msg, "\n", file = log_file, append = TRUE)
  }
  
  obj_log(paste("Starting objective function evaluation at", format(start_time)))
  
  # Determine dimensions
  Q_X <- ncol(X)
  Q_D <- dim(D)[2]
  Q_W <- ncol(W)
  Q_A <- dim(A_big)[2]
  Q_Z <- dim(Z_big)[2]
  N   <- nrow(X)
  J   <- dim(A)
  R_draws <- dim(sims_theta)[2]
  
  obj_log(sprintf("Processing data with %d observations, %d simulation draws, %d blocks", 
                  N, R_draws, Nblocks))
  
  # Set up simulation blocks
  inc <- ceiling(N / Nblocks)
  
  # Preallocate lists for likelihood and gradient for each block
  L_list <- vector("list", Nblocks)
  G_list <- vector("list", Nblocks)
  
  # Check if we should use parallel processing
  use_parallel <- !is.null(cl) & N_block > 1
  
  obj_log(sprintf("Using %s processing with %d blocks", 
                  ifelse(use_parallel, "parallel", "sequential"), Nblocks))
  
  if (use_parallel) {
    # Make variables available to worker nodes
    vars <- c("likelihood","lambda","homog")
    clusterExport(cl,  c("likelihood") , envir = environment() )
    
    # Process blocks in parallel
    obj_log("Starting parallel block processing...")
    block_results <- parLapply(cl, 1:Nblocks, function(n, param) {
      
      
      # Worker logging function
      worker_log <- function(worker_id, msg) {
        timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
        formatted_msg <- paste0("[Worker ", worker_id, " - ", timestamp, "] ", msg)
        # We use cat with file argument to avoid locking issues
        cat(formatted_msg, "\n", file = log_file, append = TRUE)
      }
      
      worker_log(n, "Starting block processing")
      
      n_min <- 1 + inc * (n - 1)
      n_max <- min(inc * n, N)
      
      
      worker_log(n, sprintf("Processing observations %d to %d (%.1f%% of data)",
                            n_min, n_max, 100 * (n_max - n_min + 1) / N))
      
      # Subset the data for block n
      A_n       <- A[n_min:n_max, ]
      A_dum_n   <- A_dum[n_min:n_max, ]
      Z_n       <- Z[n_min:n_max, ]
      E_n       <- E[n_min:n_max, ]
      X_n       <- X[n_min:n_max, ]
      D_n       <- D[n_min:n_max, , ]
      W_n       <- W[n_min:n_max, ]
      Z_big_n   <- Z_big[n_min:n_max, , ]
      A_big_n   <- A_big[n_min:n_max, , ]
      pi_i_n    <- pi_i[n_min:n_max, ]
      sims_theta_n <- sims_theta[n_min:n_max, ]
      sims_c_n     <- sims_c[n_min:n_max,  ]
      # Print each dataset to make sure we are fine 
      print(A_n)
      print(A_dum_n)
      print(Z_n)
      print(E_n)
      print(X_n)
      print(D_n)
      print(W_n)
      print(Z_big_n)
      print(A_big_n)
      print(pi_i_n)
      # Collect block data into a list
      data_n <- list(A = A_n, A_dum = A_dum_n, Z = Z_n, E = E_n, X = X_n, D = D_n, W = W_n, 
                     A_big = A_big_n, Z_big = Z_big_n, pi_i = pi_i_n, 
                     lambda = lambda, homog = homog, 
                     sims_theta = sims_theta_n, sims_c = sims_c_n,
                     log_file = log_file, worker_id = n)
      
      # Compute likelihood and gradient for the block
      res <- likelihood(param, data_n, grad = grad, worker_id = n, log_file = log_file)
      worker_log(n, "Block processing completed")
      
      return(list(L = res$L, G = res$G))
    },  omega)
    
    # Extract results from parallel blocks
    for (n in 1:Nblocks) {
      L_list[[n]] <- block_results[[n]]$L
      G_list[[n]] <- block_results[[n]]$G
    }
  } else {
    
    # Process blocks sequentially
    obj_log("Starting sequential block processing...")
    for (n in 1:Nblocks) {
      n_min <- 1 + inc * (n - 1)
      n_max <- min(inc * n, N)
      
      # Progress information
      obj_log(sprintf("Block %d/%d processing observations %d to %d (%.1f%% of data)", 
                      n, Nblocks, n_min, n_max, 100 * (n_max - n_min + 1) / N))
      
      # Subset the data for block n
      A_n       <- A[n_min:n_max, ]
      A_dum_n   <- A_dum[n_min:n_max, ]
      Z_n       <- Z[n_min:n_max, ]
      E_n       <- E[n_min:n_max, ]
      X_n       <- X[n_min:n_max, ]
      D_n       <- D[n_min:n_max, , ]
      W_n       <- W[n_min:n_max, ]
      Z_big_n   <- Z_big[n_min:n_max, , ]
      A_big_n   <- A_big[n_min:n_max, , ]
      pi_i_n    <- pi_i[n_min:n_max, ]
      sims_theta_n <- sims_theta[n_min:n_max, ]
      sims_c_n     <- sims_c[n_min:n_max,  ]
      
      # Collect block data into a list
      data_n <- list(A = A_n, A_dum = A_dum_n, Z = Z_n, E = E_n, X = X_n, D = D_n, W = W_n, 
                     A_big = A_big_n, Z_big = Z_big_n, pi_i = pi_i_n, 
                     lambda = lambda, homog = homog, 
                     sims_theta = sims_theta_n, sims_c = sims_c_n,
                     log_file = log_file, worker_id = n)
      
      # Compute likelihood and gradient for the block
      res <- likelihood(omega, data_n, grad = grad, worker_id = n, log_file = log_file)
      L_list[[n]] <- res$L
      G_list[[n]] <- res$G
      
      obj_log(sprintf("Block %d/%d completed", n, Nblocks))
    }
  }
  
  # Combine results from all blocks
  obj_log("Combining results from all blocks...")
  L_i <- do.call(c, L_list)
  # Handle any NaN or Inf values in likelihood (add small constant if needed)
  L_i[!is.finite(L_i)] <- 1e-06
  L_i[L_i <= 0] <- 1e-06
  # Construct final objective value
  Q_val <- -sum(log(L_i), na.rm = TRUE)
  obj_log(sprintf("omega: %s", paste(omega, collapse = ", ")))
  obj_log(sprintf("Objective function value: %.8f", Q_val)) #

  G_i <- NULL 
  G_val <- NULL
  # Calculate gradient if requested
  if(grad==TRUE){
    G_i <- do.call(rbind, G_list)
    G_val <- -colSums(G_i, na.rm=TRUE)
    grad_norm <- sqrt( (1/nrow(G_i)) *sum(G_val^2, na.rm=TRUE)) # norm of the average gradient is more appropriate 
    print(paste("Printing gradient:"))
    print(G_val)
    print(paste("Gradient Norm:", grad_norm))
    obj_log(sprintf("Gradient norm: %.8f", grad_norm))
  }
  end_time <- Sys.time()
  elapsed <- end_time - start_time
  obj_log(sprintf("Total elapsed time for objective function: %.2f minutes", 
                  as.numeric(elapsed)/60))
  return(list(Q = Q_val, G = G_val, G_i = G_i))
}