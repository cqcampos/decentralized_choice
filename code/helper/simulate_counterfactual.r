simulate_counterfactual <- function(dir, endyear, low_score_school, mean_school_scores, 
                                    phbao, asian, white, 
                                    yearfe, blockfe, outcome_fes, 
                                    alpha_j_math, alpha_0_X_math, alpha_m_X_math, alpha_0_math, alpha_m_math,
                                    alpha_j_ela, alpha_0_X_ela, alpha_m_X_ela, alpha_0_ela, alpha_m_ela,
                                    delta_j, delta_X, psi_D, c_X,  sd_eta,
                                    sd_theta_1, sd_theta_2, sd_theta_3,
                                    mu_theta_1, mu_theta_2, mu_theta_3,
                                    p_type_1, p_type_2, p_type_3,
                                    params, exog, round, use_torch_for_eq = FALSE, use_mean_pi = FALSE) {


  # 1. Unpack exogenous data
  X <- exog$X            # N×Q_X
  D <- exog$D            # N×Q_D×J
  W <- exog$W            # N×Q_W
  E <- exog$E
  eta <- exog$eta
  maxapps <- exog$maxapps
  last_yr <- exog$last_yr
  
  seats <- exog$seats
  seats_eq <- exog$seats_eq
  stu_type <- exog$stu_type # used for eq pi 
  param_init <- exog$param_init # replace with the last iteration values later
  
  
  
  N <- nrow(D)
  J <- dim(D)[3]
  Q_X <- ncol(X) 
  Q_D <- dim(D)[2] 
  Q_W <- ncol(W)

  # 2. Unpack some details 
  Amat            <- params$Amat  # NA×J
  pi_i            <- params$pi_i # equilibrium 
  K               <- params$K  # number of types
  homog           <- params$homog  # homogeneous flag
  preference_model<- params$preference_model
  Q_A <- ncol(Amat) 
  Q_Z <- J + 1 
  
  # Simulate shocks for nsims 
  app_shock <- -log(-log(matrix(runif(N*(J+1)), N, J+1))) 
  # post-lottery shock 
  xi <- -log(-log(matrix(runif(N*(J+1)), N, J+1))) 
  # cost shock 
  C_j <- matrix(rnorm(N*J, 0, sd_eta), N, J) 
  # Simulate latent types and theta_i 
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


  # Shift preferences for info policies 
  if(preference_model==2 | preference_model==8 | preference_model==9){ # sort of baseline info provision option
    shift_factor <- ifelse(type_draw == 1, sd_theta_1, ifelse(type_draw == 2, sd_theta_2, sd_theta_3))
    improve_draw <- runif(N)    # one Uniform(0,1) draw per obs
    theta_choice[improve_draw < 0.5] <- theta_choice[improve_draw < 0.5] +  shift_factor[improve_draw < 0.5]
  }else if(preference_model==3){
    shift_factor <- ifelse(type_draw == 1, sd_theta_1, ifelse(type_draw == 2, sd_theta_2, sd_theta_3))
    improve_draw <- (runif(N) <= 1)*(low_score_school==1)    # shock students in low-performing schools
    theta_choice[improve_draw==1] <- theta_choice[improve_draw==1] + shift_factor[improve_draw== 1]
  }
  
  
  # Calculate eq. admission probs for counterfactual policies (no need to do if considering a policy with DA )
  if( preference_model==1 | preference_model==2 | preference_model==3 | preference_model==5){
    
    if (!use_mean_pi){
      hot_start_path <- paste0(dir, "estimates/", eta, "_K", K,  "_init_param_preference_model_", preference_model,"_", last_yr, ".csv") 
      if (file.exists(hot_start_path)){
        print("Using previously optimized best-response as a starting value")
        param_init <- read.csv(hot_start_path)
        param_init <- param_init$x
      }
      
      pi_br <- calculate_eq(param = param_init,
                            eta_spec = eta,
                              seats_eq=seats_eq, pi_i=pi_i, 
                              delta_j=delta_j, delta_X=delta_X, psi_D=psi_D, c_X=c_X,
                              theta=theta_choice, C_j=C_j, 
                              X=X, D=D, W=W, A=A, E=E, 
                              stu_type=stu_type, round = round,
                              preference_model=preference_model,
                              use_torch_for_eq = use_torch_for_eq)
      
      print("Best-Response pi")
      print(pi_br)
    
    } else {
      pi_br <- read.csv(file.path(dir, "estimates", "pi_br_average.csv"))
      scenario_str <- ifelse(preference_model == 1, "Shutting Down Application Costs",
                         ifelse(preference_model == 2, "Info Provision", 
                                ifelse(preference_model == 3, "Targeted Info Provision",
                                       ifelse(preference_model == 5, "Eliminate Travel Costs",
                                              "Error: preference_model should be 1, 2, 3, or 5"))))
      
      print("Using mean(pi best response)")
      pi_br <- pi_br %>%
        filter(scenario == scenario_str) %>%
        filter(K==K) %>%
        filter(model == eta)
      
      pi_br <- pi_br$mean_pi_br
      
    }

    seats_eq$new_pi <- pi_br
    print(seats_eq, n = -1)
    # Create endyear x phbao admission probs where columns correspond to schools 
    pi_mat <- seats_eq %>% 
      ungroup() %>% 
      select(endyear, phbao, schoolcode, new_pi) %>%
      pivot_wider(
        names_from   = schoolcode,
        values_from  = new_pi,
        names_prefix = "pi_",
        values_fill  = list(new_pi = 0)   # or NA if you’d rather
      )
    # Update pi to be used in subsequent calculations 
    pi_i <- stu_type %>%
      left_join(pi_mat, by = c("endyear", "phbao")) %>% arrange(endyear, studentpseudoid)
    pi_i <- pi_i[, 4:ncol(pi_i)] 
  }
  
  # 5. Compute utilities of attending each school
  # Construct mean utilities (U_bar) 
  # We have a constant psi_j (group_matrix) so let's make a matrix of that constant 
  # delta_j <- mean(delta_j) * rep(1, J)
  delta_j_mat <- matrix(delta_j, nrow = N, ncol=length(delta_j), byrow = TRUE)
  delta_X_mat <- X %*% delta_X  
  # let's add them together and that's our U_bar 
  U_bar <- sweep(delta_j_mat, 1, delta_X_mat, FUN = "+")
  
  # CC Add distance component if we have travel costs 
  if(preference_model!=5 & preference_model!=6  & preference_model!=9){
    # Add the distance (D) component (summing over Q_D for each school j)
    D_component <- matrix(0, nrow = N, ncol = J)
    # Takes the sum of each row of D[,,j]* psi_D which is just d_{ij}*X_i*psi_D
    # It is the multiplication by X_i that initially has D[,,,] have the extra 
    # dimension 
    for (j in 1:J) {
      D_component[, j] <- rowSums(D[,,j] * matrix(psi_D, nrow = N, ncol = Q_D, byrow = TRUE))
    }
    U_bar <- U_bar + D_component
  }
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
  
  # CC For preference_model ==0,1 2, 3,5, assignments are decentralized
  if (preference_model == 0 | preference_model==1 | preference_model ==2 | preference_model==3 | preference_model==5){
    # 7. Compute application costs C
    c_X <- matrix(c_X, ncol=1, nrow=length(c_X))
    C_f_bar <- matrix( rep(exp(W %*% c_X), Q_A), ncol = Q_A) 
    if(eta =="eta_in"){
      C <- C_f_bar * exp(C_j) 
    }else{
      C <- C_f_bar + C_j 
    }
    
    # Set application costs to zero if preference_model==1 
    if (preference_model==1 ){ 
      C <-  matrix(0, nrow = N, ncol = Q_A)
    }
  
    
    # 8. Compute expected utility of each portfolio
    V <- matrix(0, nrow = N, ncol = Q_A+1)
    # Portfolio 1: Not applying anywhere
    # Expected utility is just the utility of the outside option
    V[, 1] <- Emax[, 1]/.05 + app_shock[,1] # No applications → only get the "no offers" outcome
    for (j in 1:J) {
      V[, j+1] <- (pi_i[, j] * Emax[, j+1] + (1 - pi_i[, j]) * Emax[, 1] - C[, j])/.05 + app_shock[, j+1]
    }
    
    Achoice <- max.col(V)
    # print(table(Achoice))
    A <- matrix(0, nrow = nrow(V), ncol = ncol(V))
    
    # fill a 1 in row i, column Achoice[i]
    A[cbind(seq_len(nrow(V)), Achoice)] <- 1
    A_choice <- A[, 2:ncol(A)]

    
    A_any <- Achoice >1
    
    
    # 10. Lottery: offers only to applicants
    if (FALSE){
      # 10a. Need to update lottery chances here
      gradecode <- rep(6, nrow(X))
      #Combine the ID variables with the 44 application dummies 
      A_augment <- bind_cols(
        tibble(                     # ID variables stay as ordinary columns
          endyear   = endyear,
          gradecode = gradecode,
          phbao     = phbao
        ),
        as_tibble(                  # turn the 44-column matrix into columns
          A_choice,
          .name_repair = ~ paste0("school_", seq_along(.x))   # school_1 … school_40
        )
      )
      
      ## 2.  Collapse to the desired level and add up applicants
  
      A_apps <- A_augment %>% 
        group_by(endyear, gradecode, phbao) %>%     # drop gradecode here if you
        # only want year–phbao
        summarise(
          across(                                   # choose the 40 school cols
            starts_with("school_"),                 # or use matches("school_")
            ~ sum(.x, na.rm = TRUE)
          ),
          .groups = "drop"
        )
      
      # Divide number
      seats_mat <- as.matrix(seats)[ , 4:(J+3)]      # adjust the range if you renamed
      apps_mat  <- as.matrix(A_apps)[ , 4:(J+3)]     #   or reordered columns
      school_nm <- colnames(seats_mat) 
      
      
      p_new <- seats_mat/apps_mat
      p_new[p_new > 1] <- 1 # If more seats than applicants, set p ==1 
      # If zero applicants, then 
      ##  if seats > 0, set p= 1
      ##  if seats  = 0, set p= 0   (will also catch NaN/Inf produced above) 
      zero_app <- apps_mat  == 0
      p_new[ zero_app & (seats_mat > 0) ] <- 1   # seats but nobody applied
      p_new[ seats_mat == 0 ]     <- 0 
      p_new <- cbind(as.matrix(seats)[, 1:3], p_new)
      
      
      p_new <- as.data.frame(p_new)
      stus <- cbind(endyear, phbao)
      stus <- as.data.frame(stus)
      pi_i_new <- stus %>%
        left_join(p_new, by = c("endyear", "phbao")) 
      pi_i_new <- pi_i_new[,4:ncol(pi_i_new)]  # drop ID variables
      
      #print("Pi i New")
      #print(pi_i_new)
      
      #print("A choice")
      #print(A_choice)
      Z <- (matrix(runif(N*J), N, J) < pi_i_new ) * A_choice
    }
    
    Z <- (matrix(runif(N*J), N, J) < pi_i) * A_choice
    any_offer <- rowSums(Z) > 0
    #print(mean(any_offer[Achoice!=1]))
    
    
    # 13. School selection: conditional on offers Z
    Umat     <- cbind(0,U_tilde) + (xi[,1:(J+1)]-xi[,1]) # post-lottery shock xi
    # print(rowMeans(Umat[,2:41]>0)) share of schools better than outside option during the enrollment stage
    # mean(rowMeans(Umat[,2:41]>0))
    
    choice_set <- cbind(1, Z)      # outside always in cset; only offered schools allowed
    Umat_cset <- Umat
    Umat_cset[!choice_set] <- -Inf        # Set u of schools not in cset to -inf

    Echoice_sim <- max.col(Umat_cset, ties.method = "first")

    E_new        <- diag(J+1)[Echoice_sim,]
    
    
    
  } else if (preference_model==7) {
    # This is counterfactual where students are optimally sorted based on match quality
    # Implementation of maxAllocation from outcome_analysis.do
    # Step 1: Calculate potential outcomes for each student-school pair
    # female black hispanic white poverty el eng_home lag_ela lag_math peer_q asian missing_ela missing_math
    X_outcome <- cbind(X[,1:9], X[,12], asian, X_Y[,1:2]) 
    
    X_outcome_inter <- cbind(X_outcome[, 1:11], X[, 10], X_Y[,3:7])  # Add median income te hetero and region hetero
    
    
    res_ability_math0 <- rnorm(N, 0, 1)
    Y_0_math <- as.vector(X_outcome %*% (as.matrix(alpha_0_X_math))  + theta * alpha_0_math + 
                            yearfe + res_ability_math0)
    # Y_1_j for each school j (N x J matrix)
    Y_1_j_math <-  as.vector( X_outcome_inter %*% (as.matrix(alpha_m_X_math) ) +
                              theta * (alpha_m_math) + Y_0_math) +
      matrix(rep(alpha_j_math, times = N), nrow = N, byrow = TRUE) 
    
    # Match quality for each student-school pair (N x J matrix)
    match_quality <- Y_1_j_math - as.vector(Y_0_math)
    
    # Step 2: Create a mapping from school index (1:J) to schoolcode
    # Extract unique school codes from seats_eq in order
    school_codes <- unique(seats_eq$schoolcode)
    
    # Sort schools by VAM (alpha_j) - highest quality first
    school_quality_df <- data.frame(
      school_idx = 1:J,
      school_code = school_codes[1:J],  # Assuming ordered correspondence
      quality = alpha_j_math
    )
    school_quality_df <- school_quality_df %>% arrange(desc(quality))
    
    # Step 3: Prepare seat allocation tracking
    # Initialize assignment matrix (N x J) - which school each student is assigned to
    assignment <- matrix(0, nrow = N, ncol = J)
    assigned_students <- rep(FALSE, N)
    
    # Step 4: Loop through each year
    for (yr in unique(endyear)) {
      # Students in this year
      year_mask <- (endyear == yr)
      
      # For each school in quality order
      for (school_rank in 1:nrow(school_quality_df)) {
        j <- school_quality_df$school_idx[school_rank]
        sch_code <- school_quality_df$school_code[school_rank]
        
        # Get seats available for this school-year-phbao combination
        seats_phbao1 <- seats_eq %>%
          filter(schoolcode == sch_code, endyear == yr, phbao == 1) %>%
          pull(seats) %>% 
          {if(length(.) > 0) .[1] else 0}
        
        seats_phbao0 <- seats_eq %>%
          filter(schoolcode == sch_code, endyear == yr, phbao == 0) %>%
          pull(seats) %>%
          {if(length(.) > 0) .[1] else 0}
        
        # Allocate PHBAO students
        if (seats_phbao1 > 0) {
          # Eligible students: PHBAO, in this year, not yet assigned
          eligible <- year_mask & (phbao == 1) & !assigned_students
          
          if (sum(eligible) > 0) {
            # Sort eligible students by match quality with this school
            eligible_idx <- which(eligible)
            match_scores <- match_quality[eligible_idx, j]
            top_students <- eligible_idx[order(match_scores, decreasing = TRUE)]
            
            # Allocate up to seats_phbao1 students
            n_allocate <- min(seats_phbao1, length(top_students))
            if (n_allocate > 0) {
              allocated_idx <- top_students[1:n_allocate]
              assignment[allocated_idx, j] <- 1
              assigned_students[allocated_idx] <- TRUE
            }
          }
        }
        
        # Allocate non-PHBAO students
        if (seats_phbao0 > 0) {
          # Eligible students: non-PHBAO, in this year, not yet assigned
          eligible <- year_mask & (phbao == 0) & !assigned_students
          
          if (sum(eligible) > 0) {
            # Sort eligible students by match quality with this school
            eligible_idx <- which(eligible)
            match_scores <- match_quality[eligible_idx, j]
            top_students <- eligible_idx[order(match_scores, decreasing = TRUE)]
            
            # Allocate up to seats_phbao0 students
            n_allocate <- min(seats_phbao0, length(top_students))
            if (n_allocate > 0) {
              allocated_idx <- top_students[1:n_allocate]
              assignment[allocated_idx, j] <- 1
              assigned_students[allocated_idx] <- TRUE
            }
          }
        }
      }
    }
    
    # Step 5: Create enrollment matrix E_new
    # Students assigned to a school get that column, others get outside option (column 1)
    E_new <- matrix(0, nrow = N, ncol = J + 1)
    for (i in 1:N) {
      if (any(assignment[i, ] == 1)) {
        # Student assigned to a school
        j_assigned <- which(assignment[i, ] == 1)[1]  # Get first (should be only one)
        E_new[i, j_assigned + 1] <- 1  # +1 because column 1 is outside option
      } else {
        # Student not assigned to any school - outside option
        E_new[i, 1] <- 1
      }
    }

  } else{
    # No post-application shock. Students simply rank options with "enrollment" utility greater than the outside option
    # There's no application cost, thanks to centralized admission system
    # versions 4, 6, 8, 9 should go here.
    V <- matrix(0, nrow = N, ncol = Q_A+1)
    V[,1] <- app_shock[, 1]
    V[,2:(J+1)] <- U_tilde + app_shock[, 2:(J+1)] 
    
    
    # Path where intermediary data for running DA are stored
    py_inputdir <- file.path(dir, "data/intermediate/counterfactual_matching")
    py_path <- file.path(code_dir, "code/helper/DA.py")
    
    # Run DA
    matched <- assign_da(V, stu_type, seats_eq, py_inputdir, py_path, round)
    
    # Check if the matched program gives the highest utility among ex-post feasible schools
    rol <- matched$rol 
    score <- matched$scores
    match_id <- matched$matched
    capacity <- matched$capacity
    
    ## merge score with matched, calaculate min score among matched_id 
    
    # matches + scores
    cutoffs_matches <- match_id %>%
      left_join(score, join_by(student_id == student_id))
    
    # start from capacity to keep ALL option_ids
    cutoffs <- capacity %>%
      left_join(cutoffs_matches, join_by(option_id == matched_sch_type)) %>%
      group_by(option_id) %>%
      summarise(
        capacity  = first(capacity),
        n_matched = sum(!is.na(student_id)),                        # count actual matches
        min_score = ifelse(all(is.na(score)), NA_real_,
                           min(score, na.rm = TRUE)),
        cutoff = case_when(
          capacity == 0              ~  9999.0,   # no seats
          n_matched < capacity       ~ -9999.0,   # undersubscribed
          TRUE                       ~  min_score
        ),
        .groups = "drop"
      )
    
    rol_with_cutoffs <- rol %>%
      left_join(cutoffs, join_by(option_id == option_id)) %>%
      left_join(score, join_by(id == student_id))%>%
      mutate(feas = score >= cutoff) %>%
      group_by(id) %>%
      filter(feas == 1) %>%
      mutate(u_ex_post_max = max(utility))   %>%
      mutate(ex_post_max = utility == u_ex_post_max) %>%
      left_join(match_id, join_by(id == student_id)) %>%
      mutate(matched = option_id == matched_sch_type)
    
    
    table(rol_with_cutoffs$matched, rol_with_cutoffs$ex_post_max)
  
    #last_ranked <- rol %>%
    #  group_by(id) %>%
    #  mutate(last_rank = max(rank)) %>%
    #  filter(rank == last_rank)
    
    #table(last_ranked$option_id_idx)
    

    
    E_new <- diag(J+1)[matched$matched$matched_school_id,]
    Z <- E_new[,2:(J+1)]
    
    
    

    
  }
  
  # 13. Count how many schools (endyear x costcenter x phbao) had enrollments >= seats 
  #     and applications >= seats for preference model < 4. 
  
  # Make student-type indicator matrix
  make_stu_type_ind <- function(stu_type){
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
    return(stu_type_ind_mat)
  }
  stu_type_mat <- make_stu_type_ind(stu_type)
  
  # Make matrices of seats. 1st axis is type, 2nd axis is school
  seats_mat <- matrix(seats_eq$seats, nrow = ncol(stu_type_mat))
  
  # Count enrollments for each school
  E_k <- list()
  for (k in (1:ncol(stu_type_mat))){
    E_k[[k]] <- (colSums(E_new[stu_type_mat[,k]==1, 2:(J+1)]))
  }
  E_k <- do.call(rbind, E_k)

  # Enrollments >= seats
  print("Seats")
  print(seats_mat)
  print("Enrollments")
  print(E_k)
  n_full_enroll_j <- sum(seats_mat<=E_k & seats_mat>0) 
  n_80_enroll_j  <- sum( (seats_mat*.8) <= E_k& seats_mat >0)
  
  ## CC edit 
  if (preference_model==0 | preference_model==2 | preference_model==3 | preference_model==5 ){
    # Applications >= seats
    A_k <- list()
    for (k in (1:ncol(stu_type_mat))){
      A_k[[k]] <- colSums(A_choice[stu_type_mat[,k]==1,])
    }
    A_k <- do.call(rbind, A_k)
    
    print("Actual Apps")
    print(A_k)
    n_over_sub_j <- sum(seats_mat <= A_k & seats_mat > 0 )
    
    
    print("N Offers Given")
    Z_k <- list()
    for (k in (1:ncol(stu_type_mat))){
      Z_k[[k]] <- colSums(Z[stu_type_mat[,k]==1,])
    }
    Z_k <- do.call(rbind, Z_k)
    print(Z_k)
    
    
  } else{
    n_over_sub_j <- NaN
  }
 
  # 14. Simulate outcomes Y0 and Y1
  # Y1 and Y0
  # median income not estimated because of fixed effects
  # Discrepancies in the order of variables so have to ensure this is the right arrangement 
  X_outcome <- cbind(X[,1:9], X[,12], asian, X_Y[,1:2])
  X_outcome_inter <- cbind(X_outcome[, 1:11], X[, 10], X_Y[,3:7])  # Add median income te hetero and region hetero
  
  

  res_ability_math0 <- rnorm(N, 0, 1)
  Y_0_math <- as.vector(X_outcome %*% (as.matrix(alpha_0_X_math))  + theta * alpha_0_math + 
    yearfe + res_ability_math0)
  
  Y_1_math <-  as.vector( X_outcome_inter %*% (as.matrix(alpha_m_X_math) ) +
                            theta * (alpha_m_math) + Y_0_math) +
                  matrix(rep(alpha_j_math, times = N), nrow = N, byrow = TRUE) 
  Y_1_math <- rowMeans(Y_1_math) # average across j

  res_ability_ela0 <- rnorm(N, 0, 1)
  Y_0_ela <- as.vector(X_outcome %*% (as.matrix(alpha_0_X_ela))  + theta * alpha_0_ela + 
                          yearfe + res_ability_ela0)
  Y_1_ela <-  as.vector( X_outcome_inter %*% (as.matrix(alpha_m_X_ela) ) +
                            theta * (alpha_m_ela) + Y_0_ela) +
    matrix(rep(alpha_j_ela, times = N), nrow = N, byrow = TRUE) 
  Y_1_ela <- rowMeans(Y_1_ela) # average across j

  
  
  
  # 15. Realized outcome & welfare metrics
  E_enroll <- rowSums(E_new[,2:(J+1)] >0)
  Y_math       <- Y_0_math + E_enroll*(Y_1_math - Y_0_math)
  Y_ela        <- Y_0_ela + E_enroll*(Y_1_ela - Y_0_ela)
  beta_math    <- Y_1_math - Y_0_math
  beta_ela     <- Y_1_ela - Y_0_ela
  
  # 16. Compute rates and AMTE
  #print(head(A))
  

  Z_any      <- rowSums(Z)>0
  E_any      <- E_enroll
  
  
  UC_max     <- apply(U_tilde,1,max)
  diff       <- abs(UC_max)
  sd_diff <- sqrt(var(diff, na.rm=TRUE))
  
  # threshold for marginal set at std*0.1
  threshold  <- 0.1*sd_diff
  marginal   <- (diff < threshold)*(diff>0)

  if(preference_model==0 | preference_model==1 | preference_model==2 | preference_model ==3 | preference_model==5 | preference_model==7 ){ # Decentralized
    if(preference_model==0 | preference_model==2 | preference_model==3 | preference_model==5){
      app_rate    <- mean(A_any)
      offer_rate  <- mean(Z_any[A_any])
    } else if (preference_model==7){
      app_rate    <- NaN
      offer_rate  <- NaN
    }
    attend_rate <- mean(E_any)
    attend_rate_black <- mean(E_any[X[,2]>0])
    attend_rate_hisp  <- mean(E_any[X[,3]>0])
    attend_rate_white <- mean(E_any[X[,4]>0])
    attend_rate_asian <- mean(E_any[asian>0])
    attend_rate_pov <- mean(E_any[X[,5]>0])
    

  } else{ # Centralized
    
    matched$rol %>%
      filter(rank == 1) %>%
      mutate(first_ranked = as.integer(sub("^a(\\d+)_.*", "\\1", option_id)))
    
    # Calculate share of students who applied to choice school
    first_ranked <- rol %>%
      group_by(id) %>%
      mutate(first_rank = min(rank)) %>%
      filter(rank == first_rank)
  
    n_applied <- sum(first_ranked$option_id_idx > 1)
    print("Tabulation of first ranked program idx:")
    print(table(first_ranked$option_id_idx))
    app_rate    <- n_applied/N
    
    n_choice_offer <- sum(match_id$matched_school_id > 1)
    offer_rate  <- n_choice_offer/n_applied
    attend_rate <- n_choice_offer/N
    
    n_choice_offer_black <- sum(match_id[X[,2]>0,]$matched_school_id > 1)
    n_black       <- sum(X[,2]>0)
    attend_rate_black <- n_choice_offer_black/n_black 

    n_choice_offer_hisp <- sum(match_id[X[,3] > 0,]$matched_school_id > 1)
    n_hisp              <- sum(X[,3] > 0)
    attend_rate_hisp    <- n_choice_offer_hisp / n_hisp 
    
    n_choice_offer_white <- sum(match_id[X[,4] > 0,]$matched_school_id > 1)
    n_white              <- sum(X[,4] > 0)
    attend_rate_white    <- n_choice_offer_white / n_white 
    
    n_choice_offer_asian <- sum(match_id[asian > 0,]$matched_school_id > 1)
    n_asian              <- sum(X[asian>0] > 0)
    attend_rate_asian    <- n_choice_offer_asian / n_asian 
    
    n_choice_offer_pov <- sum(match_id[X[,5] > 0,]$matched_school_id > 1)
    n_pov              <- sum(X[,5] > 0)
    attend_rate_pov    <- n_choice_offer_pov / n_pov
  }
  
  #print(app_rate)
  #print(offer_rate)
  #print(attend_rate)
  
  sd_Y_math <- sqrt(var(Y_math, na.rm=TRUE))
  tot_math         <- mean(beta_math[E_any>=1], na.rm=TRUE)
  amte_math        <- mean(beta_math[marginal==1], na.rm=TRUE)
  
  sd_Y_ela <- sqrt(var(Y_ela, na.rm=TRUE))
  tot_ela          <- mean(beta_ela[E_any>=1], na.rm=TRUE)
  amte_ela         <- mean(beta_ela[marginal==1], na.rm=TRUE)
  
  Y_black_math <- mean(Y_math[X[,2]>0], na.rm=TRUE) 
  Y_white_math <- mean(Y_math[X[,4]>0], na.rm=TRUE)
  Y_hisp_math <- mean(Y_math[X[,3]>0], na.rm=TRUE)
  Y_pov_math <- mean(Y_math[X[,5]>0], na.rm=TRUE)
  Y_nonpov_math <- mean(Y_math[X[,5]<=0], na.rm=TRUE)
  Y_asian_math <- mean(Y_math[asian>0], na.rm=TRUE)
  
  bw_gap_math <- Y_black_math - Y_white_math
  hw_gap_math <- Y_hisp_math - Y_white_math
  pov_gap_math <- Y_pov_math - Y_nonpov_math
  
  Y_black_ela <- mean(Y_ela[X[,2]>0], na.rm=TRUE)
  Y_white_ela <- mean(Y_ela[X[,4]>0], na.rm=TRUE)
  Y_hisp_ela <- mean(Y_ela[X[,3]>0], na.rm=TRUE)
  Y_pov_ela <- mean(Y_ela[X[,5]>0], na.rm=TRUE)
  Y_nonpov_ela <- mean(Y_ela[X[,5]<=0], na.rm=TRUE)
  Y_asian_ela <- mean(Y_ela[asian>0], na.rm=TRUE)
  
  bw_gap_ela <- Y_black_ela - Y_white_ela
  hw_gap_ela <- Y_hisp_ela - Y_white_ela
  pov_gap_ela <- Y_pov_ela - Y_nonpov_ela
  
  return(list(preference_model= preference_model,
              app_rate=app_rate,
              offer_rate=offer_rate,
              attend_rate=attend_rate,
              tot_math=tot_math,
              amte_math=amte_math,
              sd_Y_math = sd_Y_math,
              bw_gap_math = bw_gap_math,
              hw_gap_math = hw_gap_math,
              pov_gap_math = pov_gap_math,
              Y_black_math = Y_black_math,
              Y_white_math = Y_white_math,
              Y_hisp_math = Y_hisp_math,
              Y_pov_math = Y_pov_math,
              Y_nonpov_math = Y_nonpov_math,
              Y_asian_math = Y_asian_math,
              tot_ela=tot_ela,
              amte_ela=amte_ela,
              sd_Y_ela = sd_Y_ela,
              bw_gap_ela = bw_gap_ela,
              hw_gap_ela = hw_gap_ela,
              pov_gap_ela = pov_gap_ela,
              Y_black_ela = Y_black_ela,
              Y_white_ela = Y_white_ela,
              Y_hisp_ela = Y_hisp_ela,
              Y_pov_ela = Y_pov_ela,
              Y_nonpov_ela = Y_nonpov_ela,
              Y_asian_ela = Y_asian_ela  ,
              attend_rate_black = attend_rate_black,
              attend_rate_hisp  = attend_rate_hisp,
              attend_rate_white = attend_rate_white,
              attend_rate_asian = attend_rate_asian,
              attend_rate_pov   = attend_rate_pov,
              n_over_sub_j      = n_over_sub_j,
              n_80_enroll_j     = n_80_enroll_j,
              n_full_enroll_j   = n_full_enroll_j
              )
         )
}



# Helper function for running simulate_counterfactual for different scenarios

run_sim_cf <- function(dir, win_os, nsims, preference_model, endyear, low_score_school, mean_school_scores, 
                       phbao, asian, white, 
                       yearfe, blockfe, outcome_fes, 
                       alpha_j_math, alpha_0_X_math, alpha_m_X_math, alpha_0_math, alpha_m_math,
                       alpha_j_ela, alpha_0_X_ela, alpha_m_X_ela, alpha_0_ela, alpha_m_ela,
                       delta_j, delta_X, psi_D, c_X,  sd_eta,
                       sd_theta_1, sd_theta_2, sd_theta_3,
                       mu_theta_1, mu_theta_2, mu_theta_3,
                       p_type_1, p_type_2, p_type_3,
                       params, exog, use_torch_for_eq = FALSE, use_mean_pi = FALSE){
  
  # Declare number of clusters used for parallelization
  if (win_os){
    k <- 1  # Win OS doesn't support n_cores>1 for mclapply
  } else{
    torch_set_num_threads(1)
    Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
    
    # How many workers to use
    k <- max(1, (parallelly::availableCores() - 1))

    options(mc.cores = k)
  }
  
  # Set parent seed for mclapply
  RNGkind("L'Ecuyer-CMRG")
  set.seed(12345) 
  
  # Shut down application costs 
  params <- list(
    Amat            = Amat,
    pi_i            = pi_i,  # update later to eqm_pi[[j]],
    K               = K,
    homog           = homog,
    preference_model= preference_model
  )
  

  scenario_str <- ifelse(preference_model == 0, "baseline",
                  ifelse(preference_model == 1, "Shutting Down Application Costs",
                  ifelse(preference_model == 2, "Info Provision",
                  ifelse(preference_model == 3, "Targeted Info Provision",
                  ifelse(preference_model == 4, "DA",
                  ifelse(preference_model == 5, "Eliminate Travel Costs",
                  ifelse(preference_model == 6, "DA + No Travel Costs",
                  ifelse(preference_model == 7, "Max Allocation",
                  ifelse(preference_model == 8, "DA + Info Provision",
                  ifelse(preference_model == 9, "DA + Info Provision + No Travel Costs", NA))))))))))
  
  ids <- seq_len(nsims)
  s <- Sys.time()
  res_list <- mclapply(
    X = ids,
    FUN = function(c) {
      # prints from workers can interleave
      message(sprintf("Simulation %d of %s", c, scenario_str))
      stats <- simulate_counterfactual(
        dir,
        endyear,
        low_score_school, 
        mean_school_scores,
        phbao, asian, white,
        yearfe, blockfe, outcome_fes, 
        alpha_j_math, alpha_0_X_math, alpha_m_X_math, alpha_0_math, alpha_m_math,
        alpha_j_ela,  alpha_0_X_ela,  alpha_m_X_ela,  alpha_0_ela,  alpha_m_ela,
        delta_j, delta_X, psi_D, c_X, 
        sd_eta, 
        sd_theta_1, sd_theta_2, sd_theta_3,
        mu_theta_1, mu_theta_2, mu_theta_3, 
        p_type_1, p_type_2, p_type_3,
        params,
        exog, c,  use_torch_for_eq, use_mean_pi
      )
      stats$iter <- c  
      
      print(stats)
      return(stats)
    },
    mc.cores      = k,
    mc.preschedule = TRUE,  # good if each sim takes similar time
    mc.set.seed    = TRUE   # give each worker an independent, reproducible RNG stream
  )
  
  print(res_list)
  print(paste("Time took for", nsims, "simulations:"))
  print(Sys.time() - s)
  
  # Combine into a single data frame
  df <- do.call(rbind, lapply(res_list, as.data.frame))
  
  mean_pi_str <- ifelse(preference_model == 0 | preference_model== 4, "", ifelse(use_mean_pi, "_used_mean_pi", "_used_eq_pi"))
  
  last_yr <- exog$last_yr
  app_suffix <- ifelse(last_yr == 2013, "_2013", "")
  
  outfile <- paste0(dir, "/estimates/", eta, "_K", K , "_simulate_counterfactual_", 
                    preference_model,  mean_pi_str, app_suffix, "_maxapp_", maxapps, ".csv")
  write.csv(df, outfile, row.names = FALSE)
  
  
}