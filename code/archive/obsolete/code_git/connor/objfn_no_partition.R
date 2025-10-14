###################### Helper Function - objfn ###############################
objfn_no_partition <- function(omega, data, grad = FALSE) {
    # Start timer
    start_time <- Sys.time()
    log_file <- data$log_file
    
    K <- data$K
    
    # Compute likelihood and gradient for the block
    if (K==1){
      res <- likelihood_fast_par(omega, data, grad = grad, log_file = log_file)
    } 
    else if (K > 1){
      res <- likelihood_K_mixture(omega, data, grad = grad, log_file = log_file)
    }
    
    L_i <- res$L
    
    # Handle any NaN or Inf values in likelihood (add small constant if needed)
    L_i[!is.finite(L_i)] <- 1e-06
    L_i[L_i <= 0] <- 1e-06
    
    # Construct final objective value
    Q_val <- -sum(log(L_i), na.rm = TRUE)
    
    # Log function
    obj_log <- function(msg) {
      timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      formatted_msg <- paste0("[objfn - ", timestamp, "] ", msg)
      # We use cat with file argument to avoid locking issues
      cat(formatted_msg, "\n", file = log_file, append = TRUE)
    }
    J <- dim(data$Z)[2]
    omega_name <- c(ifelse(rep(data_max$homog, J), "Mean Utility", paste("Mean Utility", (1:J))),
                    " U x Female", " U x Black", " U x Hispanic",
                    " U x White", " U x Poverty", " U x LEP", " U x English Home",
                    " U x Lag ELA", " U x Lag Math", " U x Income",
                    " U x In Magnet", " U x Peer Q",
                    "Distance Cost", " Distance Cost x Female", " Distance Cost x Black",
                    " Distance Cost x Hispanic", " Distance Cost x White", 
                    " Distance Cost x Poverty", " Distance Cost x LEP",
                    " Distance Cost x English Home", " Distance Cost x Lag ELA",
                    " Distance Cost x Lag Math", " Distance Cost x Income",
                    " Distance Cost x In Magnet", " Distance Cost x Peer Q",
                    " Distance Cost Squared",
                    "App. Cost", " App. Cost x Female", " App. Cost x Black",
                    " App. Cost x Hispanic", " App. Cost x White",
                    " App. Cost x Poverty", " App. Cost x LEP",
                    " App. Cost x English Home", " App. Cost x Lag ELA",
                    " App. Cost x Lag Math", " App. Cost x Income",
                    " App. Cost x In Magnet", " App. Cost x Peer Q",
                    "SD eta", paste("SD theta", 1:K),   
                    ifelse(rep(K>1, K-1), paste("Alpha", 1:(K-1)), ""),
                    ifelse(rep(K>1, K-1), paste("Mu theta", 1:(K-1)), "")
    )
    
    # omega_name <- omega_name[80:84]
    
    df <- data.frame(
      Name  = omega_name,
      Value = omega,
      stringsAsFactors = FALSE
    )
    lines <- capture.output(
      print(df, row.names = FALSE, right = FALSE)
    )
    obj_log(paste(lines, collapse = "\n"))
    #obj_log(sprintf("omega: %s", paste(omega, collapse = ", ")))
    obj_log(sprintf("Objective function value: %.8f", Q_val)) #

    G_i <- NULL 
    G_val <- NULL
    # Calculate gradient if requested
    if(grad==TRUE){
        G_i <- res$G
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