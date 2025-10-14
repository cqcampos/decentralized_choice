######################### Helper Function - get_data ###########################
# Description: This helper function creates the relevant matrices and vectors
# used in the optimization procedure 
################################################################################
get_data <- function(rawdata, J, Amat) {
  # Ensure rawdata is a numeric matrix
  rawdata <- rawdata
  #storage.mode(rawdata) <- "numeric"
  
  # Initialize column pointer
  l <- 1
  N <- nrow(rawdata)
  
  # Sequentially extract variables (columns) as numeric vectors
  yearapp   <- as.numeric(rawdata$endyear)
  female    <- as.numeric(rawdata$female) 
  black     <- as.numeric(rawdata$black) 
  white     <- as.numeric(rawdata$white)
  hisp      <- as.numeric(rawdata$hispanic)
  poverty   <- as.numeric(rawdata$poverty) 
  lep       <- as.numeric(rawdata$el)
  eng_home  <- as.numeric(rawdata$eng_home) 
  current_mag <- as.numeric(rawdata$current_in_mag)
  current_peerQ <- as.numeric(rawdata$current_peer_quality)
  
  Y_m_6     <- as.numeric(rawdata$F1math) 
  Y_e_6     <- as.numeric(rawdata$F1ela)
  Y_m_7     <- as.numeric(rawdata$F2math)
  Y_e_7     <- as.numeric(rawdata$F2ela)
  Y_m_8     <- as.numeric(rawdata$F3math)
  Y_e_8     <- as.numeric(rawdata$F3ela)
  Y_m_6_8   <- as.numeric(rawdata$avg_Fmath)
  Y_e_6_8   <- as.numeric(rawdata$avg_Fela)
  
  x_m       <- as.numeric(rawdata$lag_math)
  x_e       <- as.numeric(rawdata$lag_ela)
  x_m_miss  <- as.numeric(rawdata$missing_math)
  x_e_miss  <- as.numeric(rawdata$missing_ela) 
  x_miss <- as.numeric(rawdata$median_income)
  
  Echoice   <- as.numeric(rawdata$Dchoice)
  
  # E: next (J+1) columns, ensure result is numeric matrix
  E <- as.matrix(rawdata[ , startsWith(names(rawdata), "D_")])
  #D <- as.matrix(rawdata[, l:(l + J)]); l <- l + J + 1
  
  # A: next (J-1) columns; extract as numeric matrix
  A <- as.matrix(rawdata[ , startsWith(names(rawdata), "A_")])
  Adummies <- as.matrix(rawdata[ , startsWith(names(rawdata), "A_")])
  
  # Z: next (J-1) columns
  Z <- as.matrix(rawdata[ , startsWith(names(rawdata), "Z_")])
  
  # d: next (J-1) columns
  d <- as.matrix(rawdata[ , startsWith(names(rawdata), "d_")])
  
  lat <- as.numeric(rawdata$block_lat)
  lon <- as.numeric(rawdata$block_lon) 
  
  # Index of program student applied to, as numeric
  Achoice <- as.numeric(rawdata$Achoice) 
  
  
  # Helper function for constructing column idx of applied and enrolled school
  find_col_index <- function(mat) {
    # Find column index of maximum value in each row
    result <- max.col(mat, ties.method = "first")
    # Handle rows where all values are 0
    # max.col will return 1 for all-zero rows, so we need to fix this
    row_sums <- rowSums(mat)
    result[row_sums == 0] <- -9999  
    return(result)
  }
  
  # Column index of program student 
  choice_col <- find_col_index(A)
  E_col <- find_col_index(E)
  
  # Construct X (covariates) and Y (outcomes) as numeric matrices
  X <- cbind(female, black, hisp, white, poverty, lep, eng_home, x_e, x_m, x_miss, current_mag, current_peerQ)
  X <- matrix(as.numeric(X), nrow = N)
  Y <- cbind(Y_m_6, Y_e_6, Y_m_7, Y_e_7, Y_m_8, Y_e_8, Y_m_6_8, Y_e_6_8)
  Y <- matrix(as.numeric(Y), nrow = N)
  
  # Define cohort indicators for years 2004-2008 (these will be logical; convert to numeric if needed)
  cohort <- cbind(yearapp == 2004, yearapp == 2005, yearapp == 2006, yearapp == 2007, yearapp == 2008)
  cohort <- 1 * cohort  # Convert TRUE/FALSE to 1/0
  
  Echoice <- Echoice
  return(list(N = N, A = A, Achoice = Achoice, Adummies = Adummies, Z = Z,
              choice_col = choice_col, E_col = E_col,
              E = E, Echoice = Echoice, Y = Y,
              X = X, d = d, yearapp = yearapp, cohort = cohort,
              lat = lat, lon = lon))
}
