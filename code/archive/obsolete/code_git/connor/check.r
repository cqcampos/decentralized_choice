# Calculate de for application portfolios and dVE for enrollment gradient
gA <- matrix(0, nrow = N, ncol = Q_A)
gE <- matrix(0, nrow = N, ncol = Q_A)
# No application portfolio not affected
gA[, 1] <- 0
gE[, 1] <- 0
# For each application, calculate:
for (a in 1:J) {
  # Derivative of E(u(max) \in offer)) w.r.t \delta_j
  # But it's part of the application gradients w.r.t other parameters
  gA[, a+1] <- pi_i[, a] * P_E_Z[, a+1, 2] *X[,x] 
  
}

# Calculate gA

x<- 1
# 1. pi*P_E for all applications (not including not applying)
chosenA <- rowSums(gA[, 2:Q_A]*X[,x] * A_dum , na.rm=TRUE)


res1 <- gA[, 2:Q_A] * A_dum 


###### APPLICATION STAGE GRADIENT ###### 
# Calculate dEmax directly for simplified portfolio structure
dEmax <- matrix(0, nrow = N, ncol = J+1)
Penroll <- matrix(0, nrow = N, ncol = J+1)
# No effect on the "no offers" portfolio (outside option only)
dEmax[, 1] <- 0
# For each school offer portfolio j+1
for (j in 1:J) {
  # For portfolio j+1 (offer from school j), calculate derivative of log-sum
  # Derivative equals probability of choosing school j times the covariate value (missing pi)
  # Note that this also works for the enrollment gradient 
  dEmax[, j+1] <- P_E_Z[, j+1, 2] * X[, x]
}
# Calculate gA for application portfolios
gA <- matrix(0, nrow = N, ncol = Q_A)
# No application portfolio not affected 
gA[, 1] <- 0

# For each school application
for (a in 1:J) {
  # Expected derivative from offer outcomes
  gA[, a+1] <- pi_i[, a] * dEmax[, a+1] 
}
res2 <-gA[, 2:Q_A] * A_dum

res <- data.frame(res1, res2, res1==res2, diff=res1-res2)