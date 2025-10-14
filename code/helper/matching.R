# Function that prepares data for DA matching
prep_matching <- function(V, stu_type, py_inputdir,  round, prefix = "a", id_prefix = "s") {
  stopifnot(is.matrix(V),
            all(c("endyear", "phbao") %in% names(stu_type)),
            nrow(stu_type) == nrow(V))
  
  N <- nrow(V); J <- ncol(V)
  
  # For each student, order programs by descending utility (N x J)
  idx_top <- t(apply(V, 1L, function(x) order(x, decreasing = TRUE)))
  
  # Find the position (column) where program 1 ("a1") appears in each row
  # which(..., arr.ind=TRUE) returns (row, col) of all TRUE; here exactly one per row
  rc <- which(idx_top == 1L, arr.ind = TRUE)
  pos_a1 <- integer(N); pos_a1[rc[, "row"]] <- rc[, "col"]  # length N
  
  # Build row/col indices for just the kept portion: first pos_a1[i] columns per row i
  rows_rep <- rep.int(seq_len(N), times = pos_a1)  # e.g., 1 repeated pos_a1[1] times, then 2, ...
  cols_rep <- sequence(pos_a1)                     # 1..pos_a1[1], 1..pos_a1[2], ...
  
  # Program indices for those kept (vector)
  prog_idx <- idx_top[cbind(rows_rep, cols_rep)]
  
  # Student/type pieces
  stu_suffix <- paste0("_", stu_type$endyear, "_", stu_type$phbao)
  
  # Program indices for those kept (vector)
  prog_idx <- idx_top[cbind(rows_rep, cols_rep)]
  
  # Utility values corresponding to those programs
  util_val <- V[cbind(rows_rep, prog_idx)]
  
  # Assemble long ROL
  rol <- data.frame(
    id        = paste0(id_prefix, rows_rep),
    rank      = cols_rep,
    option_id = paste0(prefix, prog_idx, stu_suffix[rows_rep]),
    option_id_idx = prog_idx,
    utility   = util_val,               
    stringsAsFactors = FALSE
  )
  
  # Create a vector of student_id 
  sids <- paste0(id_prefix, 1:N)
  
  # Create a dataframe for capacity
  # Create rows for the outside option first 
  num_types <- length(unique(stu_type$endyear))*2
  seats_outside <- seats_eq[1:num_types, ]
  seats_outside$schoolcode <- "0"
  seats_outside$seats <- rep(N, num_types)
  seats_outside$total_seats <- rep(2*N, num_types)
  seats_outside$share_seats <- rep(0.5, num_types)
  
  seats_eq <- rbind(seats_outside, seats_eq)
  option_id <- match(seats_eq$schoolcode, unique(seats_eq$schoolcode))
  option_id <- paste0(prefix, option_id, "_", seats_eq$endyear, "_", seats_eq$phbao)
  capacity <- data.frame("option_id" = option_id, 
                         "capacity" = seats_eq$seats)
  
  # Draw a single tie-breaker for each student
  scores <- data.frame("student_id" = sids,
                       "score" = runif(N))
  
  # Export ROL, student ID, and capacity 
  write.csv(rol, file.path(py_inputdir, paste0("rol_py_", round ,".csv")), row.names=FALSE) 
  write.csv(sids, file.path(py_inputdir, paste0("students_py_", round ,".csv")), row.names=FALSE) 
  write.csv(capacity, file.path(py_inputdir, paste0("capacity_", round ,".csv")), row.names=FALSE) 
  write.csv(scores, file.path(py_inputdir, paste0("scores_py_", round ,".csv")), row.names=FALSE) 
  
  return(list(rol=rol, capacity=capacity, scores = scores))
}

# Function for exporting csv files needed for DA and running it
assign_da <- function(V, stu_type, seats_eq, py_inputdir, py_path, round){
  # V: N by J+1 matrix of indirect utilities
  # stu_type : Matrix containing application year and phbao of N applicants
  # seats_eq : data.frame() containing seats available for each endyear x phbao
  # py_path: Path of DA.py
  # py_inputdir: Path where match results from DA.py will be saved
  
  # Export ROL, student ID, and capacity 
  match_inp <- prep_matching(V, stu_type, py_inputdir, round)
  
  
  # Run DA.py
  print("Running student-proposing DA algorithm to generate school choices ...")
  source_python(py_path)
  DA(py_inputdir, round)
  
  # Load match result 
  matched <- read.csv(file.path(py_inputdir, paste0("match_" , round, ".csv")))
  
  return(list(matched=matched, rol=match_inp$rol, capacity=match_inp$capacity, scores=match_inp$scores))
  
}


# Running IA with shortlist ROL and a single tie-breaker/student
ia_shortlist <- function(rol, seats_avail, tie_breaker, outside_ids = NULL) {
  # rol: N x R matrix of program IDs (1..n_prog) by rank r
  # seats_avail: length-n_prog vector of initial capacities
  # outside_ids: optional length-N vector to fill for unmatched (else NA)
  N <- nrow(rol)
  R <- ncol(rol)
  n_prog <- length(seats_avail)
  
  if (length(tie_breaker) != N){
    stop("tie_breakers must have length nrow(rol).")
  }
  matched_r    <- integer(N)            # 0 = unmatched; else rank admitted at
  seats_remain <- as.integer(seats_avail)
  
  print_match <- FALSE
  
  for (r in seq_len(R)) {
    
    # only unmatched students at this rank
    um <- matched_r == 0
    if (!any(um)) break
    
    if (print_match){
      print("=======================")
      print(paste("Round:", r))
      print("Round matched (0 if unmatched)")
      print(matched_r)
      
      print("seats remaining")
      print(seats_remain)
      
      
      print("Unmatched students")
      print((1:N)[um])
    }
    # applicants' rank-r choices (filter students who are not matched only)
    choice_r <- rol[um, r]
    
    if (!length(choice_r)) next
    
    # count unmatched applicants per program (fixed length = n_prog)
    app_r <- tabulate(choice_r, nbins = n_prog)
    
    if (print_match){
      print("Current App")
      print(app_r)
    }
    
    # undersubscribed vs oversubscribed given remaining seats
    usub <- seats_remain >= app_r
    
    ## -------- undersubscribed: admit all unmatched applicants --------
    if (any(usub)) {
      schools_u <- which(usub & app_r > 0)
      if (length(schools_u)) {
        winners_u <- which(um & rol[, r] %in% schools_u)
        if (length(winners_u)) {
          matched_r[winners_u] <- r
          # subtract exactly what we took
          taken_u <- tabulate(rol[winners_u, r], nbins = n_prog)
          seats_remain <- seats_remain - taken_u
          um[winners_u] <- FALSE
        }
      }
    }
    
    ## -------- oversubscribed: Use tiebreakers to fill in remaining capacity --------
    osub <- seats_remain < app_r
    if (any(osub)) {
      schools_o <- which(osub & app_r > 0)
      for (j in schools_o) {
        cap <- seats_remain[j]
        if (cap <= 0) next
        cands <- which(um & rol[, r] == j)
        if (!length(cands)) next
        
        top_n_idx <- function(num, n) {
          # Returns the index positions of the top n largest values in num
          if (n > length(num)) {
            stop("n cannot be larger than the length of num")
          }
          
          idx <- order(num, decreasing = TRUE)[1:n]
          return(idx)
        }
        
        
        winners_j <- cands[top_n_idx(tie_breaker[cands], cap)]
        if (print_match){
          print(paste("Winners for school", j))
          print(winners_j)
        }
        matched_r[winners_j] <- r
        seats_remain[j] <- seats_remain[j] - cap
        um[winners_j] <- FALSE
      }
    }
  }
  
  # final assignment: program ID if matched; otherwise outside_ids or NA
  matched_school <- if (is.null(outside_ids)) {
    rep(NaN, N)
  } else {
    if (length(outside_ids) != N)
      stop("outside_ids must have length nrow(rol).")
    outside_ids
  }
  m <- matched_r != 0
  if (any(m)) matched_school[m] <- rol[cbind(which(m), matched_r[m])]
  
  list(
    matched_r      = matched_r,
    matched_school = matched_school,
    seats_remain   = seats_remain
  )
}


# runs only when script is run by itself
if (sys.nframe() == 0){
  library(reticulate)
  
  py_inputdir <- "C:/users/admin/Documents/R"
  py_path <- "Z:/decentralized_choice/code/DA.py"
  source_python(py_path)
  match <- assign_da(V, stu_type, seats_eq, py_inputdir, py_path)
}