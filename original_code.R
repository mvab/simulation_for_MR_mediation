# original code by Eleanor Sanderson

datagen <- function(snps_m, beta_em, beta_mo){
  #function to generate the data
  beta_eo <- -0.4 - (beta_em*beta_mo)
  
  n=100000
  snps_e = 100
  effs_E <- (rnorm(snps_e,0.05,0.01))
  effs_M <- (rnorm(snps_m,0.05,0.01))
  
  G_e <- make_geno(n,snps_e,0.4)
  G_m <- make_geno(n,snps_m,0.4)
  v_x1 <- rnorm(n,0,1)
  v_x2 <- rnorm(n,0,1)
  C <- rnorm(n,0,1)
  
  #second sample
  G2_e <- make_geno(n,snps_e,0.4)
  G2_m <- make_geno(n,snps_m,0.4)
  v_x12 <- rnorm(n,0,1)
  v_x22 <- rnorm(n,0,1)
  C2 <- rnorm(n,0,1)
  
    Exposure <- G_e%*%effs_E +  0.8*C  + v_x1
    mediator <- G_m%*%effs_M + C + beta_em*Exposure + v_x2
    
    Exposure2 <- G2_e%*%effs_E +  0.8*C2  + v_x12
    mediator2 <- G2_m%*%effs_M + C2 + beta_em*Exposure2 + v_x22
    
  outcome <- beta_eo*Exposure2 + beta_mo*mediator2 + 0.5*C2 + rnorm(n,0,1)  
  
  #generate the summary stats
  
  G <- cbind(G_e, G_m)
  G2 <- cbind(G2_e, G2_m)
  
  MR.dat <- data.frame()
  
  for(i in 1:(snps_e+snps_m)){
    a <- summary(lm(Exposure ~ G[,i]))
    b <- summary(lm(mediator~G[,i]))
    c <- summary(lm(outcome~G2[,i]))
    MR.dat[i,"Ex_beta"] <- a$coefficient[2,1]
    MR.dat[i,"Ex_se"] <- a$coefficient[2,2]
    MR.dat[i,"Ex_p"] <- a$coefficient[2,4]
    
    MR.dat[i,"Med_beta"] <- b$coefficient[2,1]
    MR.dat[i,"Med_se"] <- b$coefficient[2,2]
    MR.dat[i,"Med_p"] <- b$coefficient[2,4]
    
    MR.dat[i,"Out_beta"] <- c$coefficient[2,1]
    MR.dat[i,"Out_se"] <- c$coefficient[2,2]
    MR.dat[i,"Out_p"] <- c$coefficient[2,4]
    
  }
  
  
  return(MR.dat)
}






## this code below generates the data and conducts the main IVW estimates
#requires packages simulateGP and dplyr
#devtools::install_github("explodecomputer/simulateGP")
library(simulateGP)
library(dplyr)

#this is what you need to add to to get the proportion mediated and loop over to run the simulation 
#get it all working with no loop and then 2 loops before you try to run it with lots of loops. 

#parameters to change:
#number of snps for the mediator
snps_m = 100

#effect of the exposure on the mediator
beta_em = -0.24

#effect of the mediator on the outcome
beta_mo = 0.07

#data generation
MR.dat <- datagen(snps_m, beta_em, beta_mo)

#analysis
ex.dat <- MR.dat[1:100,]
med.dat <- MR.dat[101:100+snps_m,]

#IVW estimation - exposure -> outcome 
lm(ex.dat$Out_beta~-1+ex.dat$Ex_beta, weights = 1/ex.dat$Out_se^2)

#IVW estimation - exposure -> mediator
lm(ex.dat$Med_beta~-1+ex.dat$Ex_beta, weights = 1/ex.dat$Med_se^2)

#IVW estimation - mediator -> outcome
lm(med.dat$Out_beta~-1+med.dat$Med_beta, weights = 1/med.dat$Out_se^2)

#IVW - MVMR estimation
lm(MR.dat$Out_beta~-1+MR.dat$Ex_beta+MR.dat$Med_beta, weights = 1/MR.dat$Out_se^2)
