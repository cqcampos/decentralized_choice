# simulate_lotteries.R
gumbel <- function(n) {               # EV1(0,1)
  -log(-log(runif(n)))
}

# Simulates one economy 
simulate_once <- function(sim,
                          exog,
                          params,
                          restrict_to_oversubscribed = TRUE,
                          scale_app = 0.05,         # scale for app-stage random utility (matches your code)
                          seed = 1L) {

  set.seed(seed + sim)
  # Sample with replacement
  N0  <- nrow(exog$X)
  idx <- sample.int(N0, N0, replace = TRUE)
  seats_df <- exog$seats_df

  X       <- exog$X[idx, , drop = FALSE]                 # N x Q_X
    W       <- exog$W[idx, , drop = FALSE]                 # N x Q_W
    D       <- exog$D[idx, , , drop = FALSE]               # N x Q_D x J
    endyear <- exog$endyear[idx]                           # length N
    phbao   <- exog$phbao[idx]                             # length N
    pi_i    <- exog$pi_i[idx, , drop = FALSE]              # N x J

  N <- nrow(X)
  Q_X <- ncol(X)
  J <- dim(D)[3]
  Q_D <- dim(D)[2]
  Q_W <- ncol(W)

  # utilities
  delta_j  <- params$delta_j            # length J
  delta_X  <- params$delta_X            # length Q_X
  psi_D    <- params$psi_D              # length Q_D

  # app cost
  c_X      <- params$c_X                # length Q_W
  sd_eta   <- params$sd_eta             # scalar; sd of app-cost idiosyncratic term
  eta <- params$eta   

  # types for theta
  K          <- params$K
  mu_theta   <- params$mu_theta         # length K
  sd_theta   <- params$sd_theta         # length K
  p_type     <- params$p_type           # length K; sums to 1

  # outcomes (math)
  alpha_j_math    <- params$alpha_j_math          # length J
  alpha_0_X_math  <- params$alpha_0_X_math        # length Q_X
  alpha_m_X_math  <- params$alpha_m_X_math        # length Q_X
  alpha_0_math    <- params$alpha_0_math          # scalar
  alpha_m_math    <- params$alpha_m_math          # scalar

  # outcomes (ela)
    alpha_j_ela     <- params$alpha_j_ela           # length J
    alpha_0_X_ela   <- params$alpha_0_X_ela         # length Q_X
    alpha_m_X_ela   <- params$alpha_m_X_ela         # length Q
    alpha_0_ela     <- params$alpha_0_ela           # scalar
    alpha_m_ela     <- params$alpha_m_ela           # scalar

  sd_eps0    <- params$sd_eps0 
  sd_eps1    <- params$sd_eps1 

  #Draw theta types and preference shocks
  type_draw <- sample.int(K, size = N, replace = TRUE, prob = p_type)
  theta_i   <- rnorm(N, mean = mu_theta[type_draw], sd = sd_theta[type_draw])

  #Application-stage shocks & post-lottery shocks
  app_shock <- matrix(gumbel(N*(J+1)), N, J+1)   # alt 1=“no apply”, 2..J+1=apply to j
  xi_post   <- matrix(gumbel(N*(J+1)), N, J+1)   # enrollment stage, option set varies by Z

  # Build utilities
  D_comp <- matrix(0, N, J)
  for (j in 1:J) {
    D_comp[, j] <- rowSums(D[, , j] * matrix(psi_D, nrow = N, ncol = Q_D, byrow = TRUE))
  }
  U_bar   <- matrix(rep(delta_j, each = N), nrow = N) +
             X %*% matrix(delta_X, ncol = 1) %*% matrix(1, nrow = 1, ncol = J) +
             D_comp
  U_tilde <- U_bar + matrix(theta_i, nrow = N, ncol = J)

  # Emax with EV1 errors (outside option = 0)
  # no-offer portfolio
  Emax0 <- rep(log(exp(0)) + 0.577215, N)

  # single-offer portfolios
  Emax_j <- matrix(0, N, J)
  for (j in 1:J) {
    Emax_j[, j] <- log(exp(0) + exp(U_tilde[, j])) + 0.577215
  }
  # Application costs
  C_bar <- as.vector(exp(W %*% matrix(c_X, ncol = 1)))           # N
  # only eta_out model
  C_j <- matrix(rnorm(N * J, 0, sd_eta), N, J)
  C    <- C_bar %o% rep(1, J) + C_j
  
  C[C < 0] <- 0  # keep non-negative

  # Application portfolio choice 
  # Value for each alternative: not applying vs applying to j
  V <- matrix(0, N, J + 1)
  V[, 1] <- Emax0/scale_app + app_shock[, 1]
  for (j in 1:J) {
    V[, j + 1] <- (pi_i[, j] * Emax_j[, j] + (1 - pi_i[, j]) * Emax0 - C[, j]) / scale_app + app_shock[, j + 1]
  }
  Achoice <- max.col(V)
  # print(table(Achoice))
  A <- matrix(0, nrow = nrow(V), ncol = ncol(V))
  
  # fill a 1 in row i, column Achoice[i]
  A[cbind(seq_len(nrow(V)), Achoice)] <- 1
  A_choice <- A[, 2:ncol(A)]
  
  
  A_any <- Achoice >1
  
  Z <- (matrix(runif(N*J), N, J) < pi_i) * A_choice
  any_offer <- rowSums(Z) > 0


  
  # 13. School selection: conditional on offers Z
  Umat     <- cbind(0,U_tilde) + (xi_post[,1:(J+1)]-xi_post[,1]) # post-lottery shock xi
  # print(rowMeans(Umat[,2:41]>0)) share of schools better than outside option during the enrollment stage
  # mean(rowMeans(Umat[,2:41]>0))
  
  choice_set <- cbind(1, Z)      # outside always in cset; only offered schools allowed
  Umat_cset <- Umat
  Umat_cset[!choice_set] <- -Inf        # Set u of schools not in cset to -inf
  
  Echoice_sim <- max.col(Umat_cset, ties.method = "first")
  
  E_new        <- diag(J+1)[Echoice_sim,]
  

  # --- Outcomes
  eps1 <- rnorm(N, 0, sd_eps1)
  eps0 <- rnorm(N, 0, sd_eps0)

  # 14. Simulate outcomes Y0 and Y1
  # Y1 and Y0
  # median income not estimated because of fixed effects
  X_outcome <- cbind(X[,1:9], X[,12])
  X_outcome <- cbind(X_outcome[,1:4], asian, X_outcome[,5:10])
  X_outcome_inter <- cbind(X[,1:10], X[,12]) # add interaction effect for median income 
  X_outcome_inter <- cbind(X_outcome_inter[,1:4], asian, X_outcome_inter[,5:11])
  
  res_ability_math <- rnorm(N, 0, 1)
  Y_1_math <-  matrix(rep(alpha_j_math, times = N), nrow = N, byrow = TRUE) +
    as.vector(X_outcome %*% (as.matrix(alpha_0_X_math)) + theta_i * (alpha_0_math +alpha_m_math) +
                X_outcome_inter %*% (as.matrix(alpha_m_X_math)) + yearfe + blockfe +
                res_ability_math)
  Y_1_math <- rowMeans(Y_1_math) # average across j
  Y_0_math <- X_outcome %*% (as.matrix(alpha_0_X_math))  + theta_i * alpha_0_math + 
    yearfe + blockfe + res_ability_math
  
  res_ability_ela <- rnorm(N, 0, 1)
  Y_1_ela <-  matrix(rep(alpha_j_ela, times = N), nrow = N, byrow = TRUE) +
    as.vector(X_outcome %*% (as.matrix(alpha_0_X_ela)) + theta_i * (alpha_0_ela + alpha_m_ela) +
                X_outcome_inter %*% (as.matrix(alpha_m_X_ela)) + yearfe_ela + blockfe_ela +
                eps1)
  Y_1_ela <- rowMeans(Y_1_ela) # average across j
  Y_0_ela <- X_outcome %*% (as.matrix(alpha_0_X_ela)) + theta_i * alpha_0_ela + 
    yearfe_ela + blockfe_ela + eps0
  

  
  # 15. Realized outcome & welfare metrics
  #print(dim(S))
  #print(dim(Y_1))
  #print(head(S))
  E_enroll <- rowSums(E_new[,2:(J+1)] >0)
  #print(S_mag)
  
  Y_math       <- Y_0_math + E_enroll*(Y_1_math - Y_0_math)
  Y_ela        <- Y_0_ela + E_enroll*(Y_1_ela - Y_0_ela)
  beta_math    <- Y_1_math - Y_0_math
  beta_ela     <- Y_1_ela - Y_0_ela
  

  # Produce long panel restricted to oversubscribed lotteries
  A_long <- as_tibble(A_choice) |>
    mutate(student = row_number(),
           endyear = endyear,
           phbao   = phbao) |>
    tidyr::pivot_longer(cols = dplyr::starts_with("V"),
                        names_to = "school", values_to = "A") |>
    mutate(school = as.integer(sub("^V", "", school))) |>
    dplyr::filter(A == 1L)

  # ---- Build school index <-> school code map from seats_df column order ----
  seat_cols    <- grep("^seats_", names(seats_df), value = TRUE)
  school_codes <- sub("^seats_", "", seat_cols)
  school_map   <- tibble::tibble(
    school = seq_along(school_codes),                 # 1..J (must match D[,,j], pi_i[,j])
    school_code = as.integer(school_codes)
  )
  stopifnot(length(seat_cols) == J)     
  
  # seats crosswalk 
  seats_key <- seats_df |>
    dplyr::select(endyear, phbao, dplyr::all_of(seat_cols)) |>
    tidyr::pivot_longer(cols = dplyr::all_of(seat_cols),
                        names_to = "seat_col", values_to = "seats") |>
    dplyr::mutate(school_code = as.integer(sub("^seats_", "", seat_col)),
                  phbao = as.integer(phbao)) |>
    dplyr::select(-seat_col) |>
    dplyr::left_join(school_map, by = "school_code") |>
    dplyr::select(endyear, phbao, school, seats, school_code)
  

  # ---- Apps, oversubscription, cell_id (by index) ----
  by_cell <- A_long |>
    dplyr::mutate(phbao = as.integer(phbao)) |>
    dplyr::count(endyear, phbao, school, name = "apps") |>
    dplyr::right_join(seats_key, by = c("endyear", "phbao", "school")) |>
    dplyr::mutate(
      apps = ifelse(is.na(apps), 0L, apps),
      oversubscribed = (apps > seats) & (seats > 0),
      cell_id = paste(endyear, phbao, school, sep = "_")
    )
  
  # ---- Offers (Z_ij) long & Enrollments (E_ij) long ----
  colnames(Z) <- NULL
  Z_long <- as_tibble(Z) |>
    dplyr::mutate(student = dplyr::row_number()) |>
    tidyr::pivot_longer(cols = dplyr::starts_with("V"),
                        names_to = "school", values_to = "Z") |>
    dplyr::mutate(school = as.integer(sub("^V", "", school)))
  
  E_long <- as_tibble(E_new[, -1]) |>
    dplyr::mutate(student = dplyr::row_number()) |>
    tidyr::pivot_longer(cols = dplyr::starts_with("V"),
                        names_to = "school", values_to = "E") |>
    dplyr::mutate(school = as.integer(sub("^V", "", school)))
  
  
  # Merge with applicants + Z + cell attributes
  # Build panel at (student, school) for applicants; attach Z, E, seats/apps/oversub
  panel <- A_long |>
    dplyr::left_join(
      by_cell |>
        dplyr::select(endyear, phbao, school, seats, apps, oversubscribed, cell_id, school_code),
      by = c("endyear", "phbao", "school")
    ) |>
    dplyr::left_join(Z_long, by = c("student", "school")) |>
    dplyr::left_join(E_long, by = c("student", "school")) |>
    dplyr::mutate(
      sim      = sim,
      applied  = A,
      offered  = Z,
      enrolled = E,
      # IMPORTANT: index Y vectors by student
      Y0_math  = Y_0_math[student],
      Y1_math  = Y_1_math[student],   # if Y_1_math is rowMeans over schools
      Y_math   = Y_math[student],
      Y0_ela   = Y_0_ela[student],
      Y1_ela   = Y_1_ela[student],
      Y_ela    = Y_ela[student],
      theta_i = theta_i[student]
    ) |>
    dplyr::select(sim, student, endyear, phbao, school, school_code,
                  cell_id, seats, apps, oversubscribed,
                  applied, offered, enrolled, 
                  Y_math, Y0_math, Y1_math, 
                  Y_ela, Y0_ela, Y1_ela, theta_i )
  
  
  if (restrict_to_oversubscribed) {
    panel <- panel |> filter(oversubscribed)
  }

  panel
}


# Puts it all together and saves to CSV file 
simulate_lottery_panel <- function(nsims,
                                   exog,
                                   params,
                                   out_file = NULL,
                                   restrict_to_oversubscribed = TRUE,
                                   scale_app = 0.05,
                                   seed = 1L) {
  dfs <- vector("list", nsims)
  for (s in 1:nsims) {
    print(paste0("Working on simulation: ", s))
    dfs[[s]] <- simulate_once(sim = s, exog = exog, params = params,
                              restrict_to_oversubscribed = restrict_to_oversubscribed,
                              scale_app = scale_app, seed = seed)
  }
  out <- bind_rows(dfs)

  if (!is.null(out_file)) {
    dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
    readr::write_csv(out, out_file)
  }
  out
}
