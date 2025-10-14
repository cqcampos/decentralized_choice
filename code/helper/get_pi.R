
########################### Helper Function - get_pi ###########################
get_pi <- function(A, Z, cohort) {
  # A, Z, and sib are matrices with dimensions (N x K)
  # cohort is a matrix with dimensions (N x ncohorts)
  N <- nrow(A)
  K <- ncol(A)
  ncohorts <- ncol(cohort)
  
  # Initialize output matrices
  pi <- matrix(0, nrow = N, ncol = K)
  pi_est <- matrix(0, nrow = K, ncol = ncohorts)
  
  # Loop over each school (k) and cohort (c)
  for (k in 1:K) {
    Atemp <- A[, k]
    Ztemp <- Z[, k]
    
    for (c in 1:ncohorts) {
      # Select individuals in the current cohort who applied (A==1) and have no sibling (sib==0)
      sel <- (Atemp == 1 & cohort[, c] == TRUE)
      if (sum(sel) > 0) {
        print(paste("This many students selected:", sum(sel)))
        Ztemp <- as.numeric(Z[which(sel==1), k])
        print(Ztemp)
        pi_est[k, c] <- mean(Ztemp, na.rm=TRUE)
        print(paste("This is admission prob:", pi_est[k, c] ) ) 
      } else {
        pi_est[k, c] <- NA
      }
      
      # For those in cohort c, assign: if sib==1 then 1; if sib==0 then pi_est
      selC <- (cohort[, c] == 1)
      if (!is.na(pi_est[k, c])) {
        pi[selC, k] <-  pi_est[k, c]
      } else {
        pi[selC, k] <- 1
      }
    }
  }
  
  return(list(pi = pi, pi_est = pi_est))
}
