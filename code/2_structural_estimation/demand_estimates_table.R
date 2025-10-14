mk_table <- function(eta, homog, K, dir){
  results <- read.csv(paste0(dir, "/estimates/in_mag_ind_", 
                            eta, 
                            "_estimates_choice_", 
                            homog, 
                            "_K", K, "_2013.csv"))
  tex_output <- file(paste0(dir, "/output/tables/demand_estimates_", 
                            eta,
                            "_",
                            homog, 
                            "_K", K, "_2013.tex"), "w")
  on.exit(close(tex_output), add = TRUE)
  ### Initiate the table ###
  writeLines (c(      "\\begin{tabular}{lccc}",
                      "\\hline \\hline",
                      "& (1) & (2) & (3) \\\\",
                      "\\hline",
                      "& & &  \\\\",
                      "& \\multicolumn{3}{c}{Panel A: Observables} \\\\",
                      "& Utility & Distance Cost & Log Cost \\\\"
                    ), con=tex_output)
  writeLines("\\cmidrule(lr){2-4}", con=tex_output)
  ### Print observables ### 
  # formatter: fixed‑format with `digits` decimals, then strip trailing 0’s
  digits <-3
  # Helper to pull out a single estimate
  pull_est <- function(prefix, cov){
    if(cov!=""){
      nm <- paste(prefix, cov)
    }else{
      nm <- prefix
    }
    val <- results$omega[results$names == nm]
    se <- results$se[results$names == nm]
    if(length(val)!=1) stop("Cannot find exactly one row for ", nm)
    val <- formatC(val, digits = digits, format = "f")
    se <- formatC(se, digits = digits, format = "f")
    list(val=val ,se=se)
  }
  # For mean utilities, report mean and sd of estimates 
  sd_mean_utilities <- sqrt( mean( 
                                  (results$omega[1:40] - mean(results$omega[1:40]) )^2 
                                    - results$se[1:40]^2 
                                  ))  
  mu_mean_utilities <- mean(results$omega[1:40])
  mu_mean_utilities_chr <- sprintf("%.3f", mu_mean_utilities)
  sd_mean_utilities_chr <- sprintf("%.3f", sd_mean_utilities)
  d <- pull_est("Distance Cost", "")
  c <- pull_est("App. Cost", "")
  writeLines( sprintf("%-12s & %6s & %6s & %6s\\\\\n",
                "Main Effects", mu_mean_utilities_chr, d[1], c[1]),
                con=tex_output)
  writeLines( sprintf(" & [%6s] & (%6s) & (%6s)\\\\\n",
                          sd_mean_utilities_chr, d[2], c[2]),
                con=tex_output)
  # Loop through rest of covariates 
  covariates <- c("Female", "Black", "Hispanic", "White",
                  "Poverty", "LEP", "English Home",
                  "Lag ELA", "Lag Math", "Income",
                  "In Magnet", "Peer Q")
  for(cov in covariates){
    u <- pull_est("U x", cov) 
    d <- pull_est("Distance Cost x", cov)
    c <- pull_est("App. Cost x", cov)
    writeLines( sprintf("%-12s & %6s & %6s & %6s\\\\\n",
                cov, u[1], d[1], c[1]),
                con=tex_output)
    writeLines(sprintf(" & (%6s) & (%6s) & (%6s)\\\\\n",
                u[2], d[2], c[2]),
                con=tex_output)
    
  }

  if(K==1){
  writeLines(c("& & &  \\\\", 
                 "& \\multicolumn{3}{c}{Panel B: Unobservables} \\\\"), 
               con=tex_output)
  theta <- pull_est("SD theta", "")
  eta <- pull_est("SD eta", "")
  writeLines("\\cmidrule(lr){2-4}", con=tex_output)
  
  # theta and eta are each a row with estimates in a 3 column multicolumn
  writeLines(paste0("$\\sigma_{\\theta}$ & \\multicolumn{3}{c}{",theta[1], "} \\\\"),
             con=tex_output)
  writeLines(paste0(" & \\multicolumn{3}{c}{(",theta[2], ")} \\\\"),
             con=tex_output)
  writeLines(paste0("$\\sigma_{\\eta}$ & \\multicolumn{3}{c}{",eta[1], "} \\\\"),
             con=tex_output)
  writeLines(paste0(" & \\multicolumn{3}{c}{(",eta[2], ")} \\\\"),
             con=tex_output)
  }
  if(K==2){
    writeLines(c("& & &  \\\\", 
                 "& \\multicolumn{3}{c}{Panel B: Unobservables} \\\\"), 
               con=tex_output)
    sd_theta1 <- pull_est("SD theta 1", "")
    sd_theta2 <- pull_est("SD theta 2", "")
    mu_theta1 <- pull_est("Mu theta 1", "")
    mu_theta2 <- pull_est("Mu theta 2", "")
    eta <- pull_est("SD eta", "")
    type1_p <- pull_est("Type 1 Probability", "")
    type2_p <- 1- as.numeric(type1_p[1])
    type2_p <- formatC(type2_p, digits = digits, format = "f")
    
    
    writeLines(paste0("& $\\mu $ & $ \\sigma$ & $Pr(K_i=k)$ \\\\"),
               con=tex_output)
    writeLines("\\cmidrule(lr){2-4}", con=tex_output)
    writeLines(paste0( "Type 1 & ", mu_theta1[1], " & ", sd_theta1[1], " & ", type1_p[1],    " \\\\ " ) ,
               con=tex_output)
    writeLines(paste0( " &  (", mu_theta1[2], ")  & (", sd_theta1[2], ") &  \\\\ " ) ,
               con=tex_output)
    
    writeLines(paste0( "Type 2 & ", mu_theta2[1], " & ", sd_theta2[1], " & ", type2_p, " \\\\ " ) ,
               con=tex_output)
    writeLines(paste0( " &  (", mu_theta2[2], ")  & (", sd_theta2[2], ") &  \\\\ " ) ,
               con=tex_output)
    
    writeLines(paste0("Cost Heterogeneity &  & ", eta[1], " &  \\\\"),
               con=tex_output)
    writeLines(paste0(" &  & (", eta[2], ") &  \\\\"),
               con=tex_output)

  }
  if(K==3){
    writeLines(c("& & &  \\\\", 
                 "& \\multicolumn{3}{c}{Panel B: Unobservables} \\\\"), 
               con=tex_output)
    sd_theta1 <- pull_est("SD theta 1", "")
    sd_theta2 <- pull_est("SD theta 2", "")
    sd_theta3 <- pull_est("SD theta 3", "")
    mu_theta1 <- pull_est("Mu theta 1", "")
    mu_theta2 <- pull_est("Mu theta 2", "")
    mu_theta3 <- pull_est("Mu theta 3", "")
    
    eta <- pull_est("SD eta", "")
    type1_p <- pull_est("Type 1 Probability", "")
    type2_p <- pull_est("Type 2 Probability", "")
    type3_p <- 1- as.numeric(type1_p[1]) - as.numeric(type2_p[1])
    type3_p <- formatC(type3_p, digits = digits, format = "f")
    
    
    writeLines(paste0("& $\\mu $ & $ \\sigma$ & $Pr(K_i=k)$ \\\\"),
               con=tex_output)
    writeLines("\\cmidrule(lr){2-4}", con=tex_output)
    writeLines(paste0( "Type 1 & ", mu_theta1[1], " & ", sd_theta1[1], " & ", type1_p[1],    " \\\\ " ) ,
               con=tex_output)
    writeLines(paste0( " &  (", mu_theta1[2], ")  & (", sd_theta1[2], ") &  \\\\ " ) ,
               con=tex_output)
    writeLines(paste0( "Type 2 & ", mu_theta2[1], " & ", sd_theta2[1], " & ", type2_p[1], " \\\\ " ) ,
               con=tex_output)
    writeLines(paste0( " &  (", mu_theta2[2], ")  & (", sd_theta2[2], ") &  \\\\ " ) ,
               con=tex_output)
    writeLines(paste0( "Type 3 & ", mu_theta3[1], " & ", sd_theta3[1], " & ", type3_p, " \\\\ " ) ,
               con=tex_output)
    writeLines(paste0( " &  (", mu_theta3[2], ")  & (", sd_theta3[2], ") &  \\\\ " ) ,
               con=tex_output)
    
    writeLines(paste0("Cost Heterogeneity &  & ", eta[1], " &  \\\\"),
               con=tex_output)
    writeLines(paste0(" &  & (", eta[2], ") &  \\\\"),
               con=tex_output)
    
  }
    # Close the table
  writeLines(c("& & &  \\\\", "\\hline"), con=tex_output)
  writeLines(c("\\hline", "\\end{tabular}"), con=tex_output)
  
}
dir <- "/Volumes/lausd/decentralized_choice/"
#mk_table(eta="eta_in", homog="hetero", K=1, dir=dir)
mk_table(eta="eta_out", homog="hetero", K=1, dir=dir)
#mk_table(eta="eta_in", homog="hetero", K=2, dir=dir)
mk_table(eta="eta_out", homog="hetero", K=2, dir=dir)
#mk_table(eta="eta_in", homog="hetero", K=3, dir=dir)
#mk_table(eta="eta_out", homog="hetero", K=3, dir=dir)

