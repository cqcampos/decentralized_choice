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
likelihood <- function(omega, data, grad = FALSE, log_file = NULL, posterior = FALSE) {
  if (is.null(log_file) && !is.null(data$log_file)) log_file <- data$log_file
  lk_log <- function(msg) {
    if (!is.null(log_file)) {
      timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      cat(paste0(" (likelihood) - ", timestamp, "] ", msg, "\n"),
          file = log_file, append = TRUE)
    }
  }
  
  ##############################################################################
  # Data
  ##############################################################################
  likelihood_i <- NULL
  G_i <- NULL
  
  A          <- data$A
  A_dum      <- data$A_dum
  Achoice    <- data$Achoice
  applied_in <- data$applied
  Z          <- data$Z
  X          <- data$X
  E          <- data$E
  W          <- data$W
  A_big      <- data$A_big
  Z_big      <- data$Z_big
  pi_i       <- data$pi_i
  lambda     <- data$lambda
  homog      <- data$homog
  sims_theta <- data$sims_theta   # N x R
  sims_c     <- data$sims_c       # N x J x R
  cl         <- data$cl
  choice_col <- data$choice_col
  E_col      <- data$E_col
  eta_out    <- data$eta_out
  n_cores    <- data$num_cores
  D          <- data$D
  K          <- data$K
  
  rm(data); gc()
  
  Q_X     <- ncol(X)
  Q_D     <- dim(D)[2]
  Q_W     <- ncol(W)
  Q_A     <- dim(A_big)[2]
  Q_Z     <- ncol(Z)
  N       <- nrow(X)
  J       <- ncol(Z)
  R_draws <- ncol(sims_theta)
  
  ##############################################################################
  # Parameters
  ##############################################################################
  l_param <- 1
  
  # 1) delta_j
  if (homog == 1) {
    delta_j <- rep(omega[l_param], J)
    l_param <- 2                      # <-- EXACTLY as in the snippet you sent
  } else {
    delta_j <- omega[l_param:(l_param + J - 1)]
    l_param <- l_param + J
  }
  
  # 2) beta_X
  beta_X <- matrix(omega[l_param:(l_param + Q_X - 1)], nrow = Q_X, ncol = 1)
  l_param <- l_param + Q_X
  
  # 3) psi_D
  psi_D <- matrix(omega[l_param:(l_param + Q_D - 1)], nrow = Q_D, ncol = 1)
  l_param <- l_param + Q_D
  
  # 4) c_f (fixed application costs)
  c_f <- omega[l_param:(l_param + Q_W - 1)]
  l_param <- l_param + Q_W
  
  # 5) sigma_c (SD of eta_i in costs) — exponentiated
  sigma_c <- exp(omega[l_param])
  l_param <- l_param + 1
  
  # 6) sigma_theta (length K) — exponentiated
  sigma_theta <- exp(omega[l_param:(l_param + K - 1)])
  l_param <- l_param + K
  
  # 7) alpha_k (K-1 free, last fixed to 0)
  if (K > 1) {
    alpha_k <- c(omega[l_param:(l_param + K - 2)], 0)
    l_param <- l_param + (K - 1)
  } else {
    alpha_k <- 0
  }
  
  #    softmax -> p_k   (use same math as the snippet; no stabilization to stay exact)
  exp_alpha_k     <- exp(alpha_k)
  sum_exp_alpha_k <- sum(exp_alpha_k)
  p_k <- exp_alpha_k / sum_exp_alpha_k    # mixture weights (length K)
  
  # 8) mu_theta (K-1 free, last implied by zero-mean restriction)
  if (K > 1) {
    mu_free  <- omega[l_param:(l_param + K - 2)]
    l_param  <- l_param + (K - 1)
    mu_theta <- c(mu_free, -sum(p_k[1:(K - 1)] * mu_free) / p_k[K])
  } else {
    mu_theta <- 0
  }
  
  if (posterior) {
    log_prior_r <- NULL
    theta_r     <- NULL
  }
  
  ##############################################################################
  # Deterministic parts
  ##############################################################################
  group_matrix <- matrix(rep(delta_j, each = N), nrow = N)
  fixed_effect <- X %*% beta_X
  U_bar <- sweep(group_matrix, 1, fixed_effect, "+")
  D_component <- matrix(0, N, J)
  for (j in 1:J) {
    D_component[, j] <- rowSums(D[,, j] * matrix(psi_D, N, Q_D, byrow = TRUE), na.rm = TRUE)
  }
  U_bar <- U_bar + D_component
  
  applied <- rowSums(A) > 0
  C_f_bar <- matrix(rep(exp(W %*% c_f), Q_A), ncol = Q_A)
  C_f_bar[, 1] <- 0
  
  hasOffers <- rowSums(Z) > 0
  enrolled  <- rowSums(E[, 2:Q_A]) > 0
  noEnrollmentStageChoice  <- (applied == 0 | (applied == 1 & hasOffers == 0))
  applied_chose_school     <- (applied == 1 & hasOffers == 1 & enrolled == 1)
  
  ##############################################################################
  # Simulated likelihood
  ##############################################################################
  Q_params <- length(omega)
  start_time <- Sys.time()
  
  Pbar_i_k <- matrix(0, N, K)
  compute_grad <- grad && (K == 1)
  
  res_store <- NULL  ### keep per-draw results for K==1 gradients
  
  l_and_grad_r_k <- function(r, k_idx, mu_k, sigma_k) {
    sims_theta_r <- matrix(sims_theta[, r], ncol = 1)
    theta <- mu_k + sigma_k * sims_theta_r
    
    eta_ij <- sigma_c * sims_c[,, r]
    
    eta_mat <- cbind(matrix(0, N, 1), eta_ij)
    if (eta_out) {
      C_total <- C_f_bar + eta_mat
    } else {
      C_total <- C_f_bar * exp(eta_mat)
    }
    C_total <- cbind(matrix(0, N, 1), C_total[, 2:ncol(C_total)])
    
    U_tilde <- U_bar + matrix(rep(theta, J), nrow = N, ncol = J)
    exp_U_0 <- 1
    exp_U_1 <- exp(U_tilde)
    
    EULER <- 0.5772156649
    Emax <- matrix(0, N, J + 1)
    Emax[, 1] <- EULER
    for (j in 1:J) {
      Emax[, j + 1] <- log(exp_U_0 + exp_U_1[, j]) + EULER
    }
    
    e <- matrix(0, N, Q_A)
    e[, 1] <- Emax[, 1] / lambda
    for (j in 1:J) {
      e[, j + 1] <- (1 / lambda) * pi_i[, j] * Emax[, j + 1] +
        (1 / lambda) * (1 - pi_i[, j]) * Emax[, 1] -
        (C_total[, j + 1]) / lambda
    }
    e <- exp(e)
    
    row_sums <- rowSums(e)
    P_A <- e / row_sums
    A_dum_full <- cbind(1 - rowSums(A_dum), A_dum)
    P_A_i <- rowSums(P_A * A_dum_full)
    
    numerator   <- Z * exp_U_1
    denominator <- 1 + rowSums(numerator)
    P_E <- matrix(0, N, J + 1)
    no_offers <- rowSums(Z) == 0
    P_E[no_offers, 1]  <- 1
    P_E[!no_offers, 1] <- 1 / denominator[!no_offers]
    for (j in 1:J) {
      P_E[!no_offers, j + 1] <- numerator[!no_offers, j] / denominator[!no_offers]
    }
    P_E_i <- rowSums(P_E * E)
    
    likelihood_r <- P_A_i * P_E_i
    
    if (!compute_grad) {
      return(list(likelihood_r = likelihood_r)) # not gradients for k>1 right now.
    }
    
    ########################## GRADIENTS (K == 1 only) #########################
    P_E_Z <- array(0, dim = c(N, J + 1, 2))
    P_E_Z[, 1, 1] <- 1
    for (j in 1:J) {
      P_E_Z[, j + 1, 1] <- exp_U_0 / (exp_U_0 + exp_U_1[, j])
      P_E_Z[, j + 1, 2] <- exp_U_1[, j] / (exp_U_0 + exp_U_1[, j])
    }
    
    P_A_short <- P_A[, 2:Q_A]
    
    gA <- matrix(0, N, Q_A)
    gE <- matrix(0, N, Q_A)
    for (a in 1:J) {
      gA[, a + 1] <- pi_i[, a] * P_E_Z[, a + 1, 2]
      gE[, a + 1] <- P_E_Z[, a + 1, 2]
    }
    chosenA <- rowSums(gA[, 2:Q_A] * A_dum)
    gA <- gA[, 2:Q_A]
    chosenE <- rowSums(gE[, 2:Q_A] * A_dum)
    
    if (homog == 1) {
      dlog_P_A_delta <- (chosenA - rowSums(P_A_short * gA)) / lambda
      dlog_P_E_delta <- (ifelse(applied_chose_school, 1, 0) - chosenE) * (1 - noEnrollmentStageChoice)
    } else {
      dlog_P_A_delta <- matrix(0, N, J)
      dlog_P_E_delta <- matrix(0, N, J)
      for (j in 1:J) {
        dlog_P_A_delta[, j] <- (1 / lambda) * chosenA * (choice_col == j) -
          (1 / lambda) * P_A_short[, j] * gA[, j]
        dlog_P_E_delta[, j] <- (applied_chose_school - chosenE) * (choice_col == j) * (1 - noEnrollmentStageChoice)
      }
    }
    
    dlog_P_A_X <- matrix(0, N, Q_X)
    dlog_P_E_X <- matrix(0, N, Q_X)
    for (x in 1:Q_X) {
      dlog_P_A_X[, x] <- (1 / lambda) * chosenA * X[, x] -
        rowSums((1 / lambda) * gA * X[, x] * P_A_short)
      dlog_P_E_X[, x] <- (ifelse(applied_chose_school, 1, 0) - chosenE) * X[, x] * (1 - noEnrollmentStageChoice)
    }
    
    dlog_P_A_D <- matrix(0, N, Q_D)
    dlog_P_E_D <- matrix(0, N, Q_D)
    for (d in 1:Q_D) {
      relCovariate <- rowSums(D[, d, ] * A_dum)
      dlog_P_A_D[, d] <- (1 / lambda) * chosenA * relCovariate -
        rowSums((1 / lambda) * gA * D[, d, ] * P_A_short)
      dlog_P_E_D[, d] <- (ifelse(applied_chose_school, 1, 0) - chosenE) * relCovariate * (1 - noEnrollmentStageChoice)
    }
    
    if (eta_out) {
      dC <- eta_ij
    } else {
      dC <- C_total[, 2:ncol(C_total)] * eta_ij
    }
    rm(eta_ij); gc()
    chosendC <- rowSums(dC * A_dum)
    dlog_P_A_sd_eta <- -((1 / lambda) * chosendC - rowSums((1 / lambda) * P_A_short * dC))
    
    # theta sd (K==1)
    dU <- theta
    dlog_P_A_theta <- (1 / lambda) * chosenA * dU -
      rowSums((1 / lambda) * P_A_short * gA * dU[, rep(1, ncol(gA))])
    dlog_P_E_theta <- (ifelse(applied_chose_school, dU, 0) - chosenE * dU) *
      (1 - noEnrollmentStageChoice)
    
    dlog_P_A_W <- matrix(0, N, Q_W)
    if (eta_out) {
      for (w in 1:Q_W) {
        dCw <- C_f_bar[, 2:Q_A] * W[, w]
        dCapplied <- rowSums(dCw * A_dum)
        dlog_P_A_W[, w] <- -((1 / lambda) * dCapplied - rowSums((1 / lambda) * P_A_short * dCw))
      }
    } else {
      for (w in 1:Q_W) {
        dCw <- C_total * W[, w]
        dCw <- dCw[, -1]
        dCapplied <- rowSums(dCw * A_dum)
        dlog_P_A_W[, w] <- -((1 / lambda) * dCapplied - rowSums((1 / lambda) * P_A_short * dCw))
      }
    }
    
    dlog_P_E_sd_eta <- matrix(0, N, 1)
    dlog_P_E_W      <- matrix(0, N, Q_W)
    
    dlog_P_A <- cbind(dlog_P_A_delta, dlog_P_A_X, dlog_P_A_D,
                      dlog_P_A_sd_eta, dlog_P_A_theta, dlog_P_A_W)
    dlog_P_E <- cbind(dlog_P_E_delta, dlog_P_E_X, dlog_P_E_D,
                      dlog_P_E_sd_eta, dlog_P_E_theta, dlog_P_E_W)
    
    dlikelihood_sim_r <- (dlog_P_A + dlog_P_E) * likelihood_r
    
    list(likelihood_r = likelihood_r,
         dlikelihood_sim_r = dlikelihood_sim_r)
  }
  
  # ------------------------------- main loops -------------------------------- #
  for (k in 1:K) {
    mu_k    <- mu_theta[k]
    sigma_k <- sigma_theta[k]
    
    res_k <- parallel::mclapply(
      1:R_draws,
      function(r) l_and_grad_r_k(r, k, mu_k, sigma_k),
      mc.cores = n_cores
    )
    
    likelihood_sim_k <- do.call(cbind, lapply(res_k, `[[`, "likelihood_r"))
    Pbar_i_k[, k] <- rowMeans(likelihood_sim_k, na.rm = TRUE)
    
    if (compute_grad) res_store <- res_k  ### save for aggregation
  }
  
  total_time <- Sys.time() - start_time
  lk_log(sprintf("Completed all %d simulation draws over %d types in %.2f minutes",
                 R_draws, K, as.numeric(total_time) / 60))
  
  # mixture likelihood
  likelihood_i <- as.vector(Pbar_i_k %*% p_k)
  
  # log the pieces you asked for
  for (kk in 1:K) {
    lk_log(sprintf("sum(Pbar_i_k[,%d]) = %.6e", kk, sum(Pbar_i_k[, kk], na.rm = TRUE)))
  }
  lk_log(sprintf("p_k = [%s]", paste(format(p_k, digits = 6), collapse = ", ")))
  lk_log(sprintf("sum(likelihood_i) = %.6e", sum(likelihood_i, na.rm = TRUE)))
  lk_log(sprintf("sum(-log(likelihood_i)) = %.6e", sum(-log(likelihood_i), na.rm = TRUE)))
  
  # diagnostics (you can keep or remove the prints)
  print(summary(likelihood_i))
  print(quantile(likelihood_i,
                 probs = c(0.01, 0.05, 0.10, 0.15, 0.20, 0.40, 0.5,
                           0.60, 0.80, 0.90, 0.95, 0.99),
                 na.rm = TRUE))
  
  threshold_count <- sum(likelihood_i < 1e-6, na.rm = TRUE)
  percent_threshold <- round(100 * threshold_count / length(likelihood_i), 2)
  lk_log(sprintf("Low likelihood warning: %d observations (%g%%) hit minimum threshold (1e-6)",
                 threshold_count, percent_threshold))
  lk_log(sprintf("Pre-calculation check: %d NaNs in likelihood_i", sum(is.na(likelihood_i))))
  lk_log(sprintf("Likelihood summary - Min: %.3e, Q1: %.3e, Median: %.3e, Mean: %.3e, Q3: %.3e, Max: %.3e",
                 min(likelihood_i, na.rm = TRUE), quantile(likelihood_i, 0.25, na.rm = TRUE),
                 median(likelihood_i, na.rm = TRUE), mean(likelihood_i, na.rm = TRUE),
                 quantile(likelihood_i, 0.75, na.rm = TRUE), max(likelihood_i, na.rm = TRUE)))
  
  # gradients
  if (grad && K == 1) {
    stopifnot(!is.null(res_store))
    
    likelihood_sim <- do.call(cbind, lapply(res_store, `[[`, "likelihood_r"))
    
    array_dim <- c(N, Q_params, R_draws)
    dlikelihood_sim <- array(
      unlist(lapply(res_store, `[[`, "dlikelihood_sim_r")),
      dim = array_dim
    )
    
    dlikelihood_i <- apply(dlikelihood_sim, c(1, 2), mean, na.rm = TRUE)
    likelihood_i_mat <- matrix(likelihood_i, nrow = N, ncol = ncol(dlikelihood_i))
    G_i <- dlikelihood_i / likelihood_i_mat
    
  } else if (grad && K > 1) {
    # implement k>1 gradients.
  }
  
  list(
    L = likelihood_i,
    G = G_i,
    posterior_mu = NULL,
    posterior_sd = NULL,
    p_k = p_k,
    mu_theta = mu_theta,
    sigma_theta = sigma_theta
  )
}
