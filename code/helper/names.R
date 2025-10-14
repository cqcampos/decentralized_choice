param_names <- function(K, last_yr){
  J <- ifelse(last_yr == 2013, 53, 40)
  if(K==1){
    mean_util <- sprintf("Mean Utility %d", 1:J)
    
    # 2) list the remaining interaction and cost terms
    others <- c(
      "U x Female",       "U x Black",        "U x Hispanic",     "U x White",
      "U x Poverty",      "U x LEP",          "U x English Home", "U x Lag ELA",
      "U x Lag Math",     "U x Income",       "U x In Magnet",    "U x Peer Q",
      
      "Distance Cost",    "Distance Cost x Female", "Distance Cost x Black",
      "Distance Cost x Hispanic",    "Distance Cost x White",
      "Distance Cost x Poverty",     "Distance Cost x LEP",
      "Distance Cost x English Home","Distance Cost x Lag ELA",
      "Distance Cost x Lag Math",    "Distance Cost x Income",
      "Distance Cost x In Magnet",   "Distance Cost x Peer Q",
      "Distance Cost Squared",
      "App. Cost",        "App. Cost x Female", "App. Cost x Black",
      "App. Cost x Hispanic","App. Cost x White","App. Cost x Poverty",
      "App. Cost x LEP",  "App. Cost x English Home", "App. Cost x Lag ELA",
      "App. Cost x Lag Math","App. Cost x Income","App. Cost x In Magnet",
      "App. Cost x Peer Q", 
      "SD eta",           "SD theta"
    )
    names <- c(mean_util, others)

  }
  if(K==2){
    names <- c(
      sprintf("Mean Utility %d", 1:J),
      
      # Interaction terms with U
      "U x Female", "U x Black", "U x Hispanic", "U x White",
      "U x Poverty", "U x LEP", "U x English Home",
      "U x Lag ELA", "U x Lag Math", "U x Income",
      "U x In Magnet", "U x Peer Q",
      
      # Distance cost and its interactions
      "Distance Cost",
      "Distance Cost x Female", "Distance Cost x Black", "Distance Cost x Hispanic",
      "Distance Cost x White", "Distance Cost x Poverty", "Distance Cost x LEP",
      "Distance Cost x English Home", "Distance Cost x Lag ELA",
      "Distance Cost x Lag Math", "Distance Cost x Income",
      "Distance Cost x In Magnet", "Distance Cost x Peer Q",
      "Distance Cost Squared",
      
      # Application cost and its interactions
      "App. Cost",
      "App. Cost x Female", "App. Cost x Black", "App. Cost x Hispanic",
      "App. Cost x White", "App. Cost x Poverty", "App. Cost x LEP",
      "App. Cost x English Home", "App. Cost x Lag ELA",
      "App. Cost x Lag Math", "App. Cost x Income",
      "App. Cost x In Magnet", "App. Cost x Peer Q",
      
      # Standard deviations and mixture parameters
      "SD eta", "SD theta 1", "SD theta 2",
      "Alpha 1", "Mu theta 1",
      "Type 1 Probability",
      "Mu theta 2"
    )
  }
  if(K==3){
    names <- c(
      sprintf("Mean Utility %d", 1:J),
      
      # Interaction terms with U
      "U x Female", "U x Black", "U x Hispanic", "U x White",
      "U x Poverty", "U x LEP", "U x English Home",
      "U x Lag ELA", "U x Lag Math", "U x Income",
      "U x In Magnet", "U x Peer Q",
      
      # Distance cost and its interactions
      "Distance Cost",
      "Distance Cost x Female", "Distance Cost x Black", "Distance Cost x Hispanic",
      "Distance Cost x White", "Distance Cost x Poverty", "Distance Cost x LEP",
      "Distance Cost x English Home", "Distance Cost x Lag ELA",
      "Distance Cost x Lag Math", "Distance Cost x Income",
      "Distance Cost x In Magnet", "Distance Cost x Peer Q",
      "Distance Cost Squared",
      
      # Application cost and its interactions
      "App. Cost",
      "App. Cost x Female", "App. Cost x Black", "App. Cost x Hispanic",
      "App. Cost x White", "App. Cost x Poverty", "App. Cost x LEP",
      "App. Cost x English Home", "App. Cost x Lag ELA",
      "App. Cost x Lag Math", "App. Cost x Income",
      "App. Cost x In Magnet", "App. Cost x Peer Q",
      
      # Standard deviations and mixture parameters
      "SD eta", "SD theta 1", "SD theta 2", "SD theta 3",
      "Alpha 1", "Alpha 2",
      "Mu theta 1", "Mu theta 2",
      "Type 1 Probability", "Type 2 Probability",
      "Mu theta 3"
    )
  }
  return(names)
}