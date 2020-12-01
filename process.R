library(parallel)
library(MASS)
library(simulateGP)
library(dplyr)
library(readr)
library(lubridate)
source("original_code.R")


if(!exists("args") || class(args)!="character"){
  args <- commandArgs(trailingOnly = TRUE)
}

if (length(args)==0) {
  
  cat("MR results simulation script
      
      Usage: Rscript <this_script>.R <numCores> <numIter> [<data_path>]
      * <numCores> : number of cores to use to run analyses in paralell
      * <numIter> : number of simulation iterations to run (1 iter takes 1 min on 1 core)
      * <data_path>: data directory; if not supplied, will assumed that data folder in the project folder should be used
      
      ")
  stop("Must provide at least 6 arguments")
  
} else if ( length(args) == 3 ){
  
  numCores <- args[1] 
  numIter <- args[2] 
  data_path <- args[3]
  
} else {
  numCores <- args[1]
  numIter <- args[2] 
  data_path <- "todo"
}

time0 <- gsub(" ", "_", now())
print(paste0("====== Time 0: ", time0))

input <-read_csv(paste0(data_path, "/simulation_input.csv"))

# just trying for first mediator
i = 1

#number of snps for the mediator
snps_m = as.numeric(input[i, "mediator_SNPs"])

#effect of the exposure on the mediator
beta_em = as.numeric(input[i, "effect_exp_med"])

#effect of the mediator on the outcome
beta_mo = as.numeric(input[i, "effect_med_out"])

### FUNCTIONS

simulate_results <- function(n_iter, snps_m, beta_em, beta_mo){
  
  #data generation
  MR.dat <- datagen(snps_m, beta_em, beta_mo)

  #analysis
  ex.dat <- MR.dat[1:100, ]
  med.dat <- MR.dat[101:as.numeric(100+snps_m), ] 
  
  #IVW estimation - exposure -> outcome 
  eo <- lm(ex.dat$Out_beta ~ -1+ex.dat$Ex_beta, weights = 1/ex.dat$Out_se^2)
  eo_beta <- summary(eo)$coefficients[1] 
  eo_se <- summary(eo)$coefficients[2] 
  
  #IVW estimation - exposure -> mediator
  em <- lm(ex.dat$Med_beta ~ -1+ex.dat$Ex_beta, weights = 1/ex.dat$Med_se^2)
  em_beta <- summary(em)$coefficients[1] 
  em_se <- summary(em)$coefficients[2] 
  
  #IVW estimation - mediator -> outcome
  mo <- lm(med.dat$Out_beta ~ -1+med.dat$Med_beta, weights = 1/med.dat$Out_se^2)
  mo_beta <- summary(mo)$coefficients[1] 
  mo_se <- summary(mo)$coefficients[2] 
  
  #IVW - MVMR estimation
  emo <- lm(MR.dat$Out_beta~-1+MR.dat$Ex_beta+MR.dat$Med_beta, weights = 1/MR.dat$Out_se^2)
  emo_betas <- summary(emo)$coefficients[,1] %>% as.numeric()
  emo_ses <- summary(emo)$coefficients[,2] %>% as.numeric()
  
  
  
  output_vector <-  data.frame(eo_beta, em_beta, mo_beta, emo_betas[1], emo_betas[2],
                               eo_se, em_se, mo_se, emo_ses[1], emo_ses[2])  
  
  output_df<-rbind(data.frame(), output_vector)
  
  print(paste("Simulation iteration number", n_iter, "completed at:",  now()  ))
  
  return(output_df)
}


difference_method <- function(EO_beta_total, EO_beta_direct){
  # calculate indirect effect beta
  
  # INDIRECT = TOTAL (of exposure, univ) - DIRECT (of exposure, mvmr)
  indirect_beta = EO_beta_total - EO_beta_direct
  return(indirect_beta)
}

product_method <- function(EM_beta, MO_beta){
  # calculate indirect effect beta
  
  # method 1
  # INDIRECT = TOTAL (exposure -> mediator) x TOTAL (mediator -> outcome)
  # method 2
  # INDIRECT = TOTAL (exposure -> mediator) x DIRECT (of mediator -> outcome , mvmr) 
  
  indirect_beta =  EM_beta * MO_beta
  return(indirect_beta)
}


delta_method <- function(EM_beta, EM_se, MO_beta, MO_se){

  # SE of INDIRECT effect (applied with product method) 
  delta_se = sqrt((MO_beta^2 * EM_se^2) + (EM_beta^2 * MO_se^2))
  return(delta_se)
}

propagation_of_errors_method <- function(total_se, direct_se){
  
  # SE of INDIRECT effect (applied with difference method) = sqrt(SE TOTAL^2 + SE DIRECT^2) 
  # SE of INDIRECT effect (applied with product method) = sqrt(SE EM^2 + SE MO^2) 
  
  indirect_se = sqrt(total_se^2 + direct_se^2)
  return(indirect_se)
  
}



### MAIN code

# simulate results and produce summary effects `numIter` times
n_iter = seq(1:numIter)
results_list <- mclapply(n_iter, simulate_results, mc.cores = numCores,
                                snps_m = snps_m,
                                beta_em = beta_em, 
                                beta_mo = beta_mo)


# convert list of vectors to df
results <- bind_rows(lapply(results_list, as.data.frame.list))
colnames(results) <- c("EO_total_beta", "EM_total_beta", "MO_total_beta", "EO_direct_beta", 'MO_direct_beta',
                       "EO_total_se", "EM_total_se", "MO_total_se", "EO_direct_se", 'MO_direct_se')


# calculate indirect beta using difference and product method (x2)
results <- results %>% 
      # mediation using difference method and PoE for SE calulation
      mutate(indirect_b_difference = difference_method(EO_total_beta, EO_direct_beta)) %>% 
      mutate(indirect_se_difference_PoE = propagation_of_errors_method(EO_total_se, EO_direct_se)) %>% 
  
      # mediation using Product method (total effect of two steps) and both PoE and Delta for SE calculation
      mutate(indirect_b_product_v1 = product_method(EM_total_beta, MO_total_beta)) %>% 
      mutate(indirect_se_product_v1_PoE = propagation_of_errors_method(EM_total_se, MO_total_se)) %>% 
      mutate(indirect_se_product_v1_delta = delta_method(EM_total_beta,EM_total_se,MO_total_beta, MO_total_se)) %>% 
                                    
      # mediation using Product method (total effect of EM and direct effect of MO) and both PoE and Delta for SE calculation
      mutate(indirect_b_product_v2 = product_method(EM_total_beta, MO_direct_beta)) %>% 
      mutate(indirect_se_product_v2_PoE = propagation_of_errors_method(EM_total_se, MO_direct_se)) %>% 
      mutate(indirect_se_product_v2_delta = delta_method(EM_total_beta,EM_total_se,MO_direct_beta, MO_direct_se))
        

# also save input parameters to file
results$param_snps_m <- snps_m
results$param_beta_em <- beta_em
results$param_beta_mo <- beta_mo


time1 <- gsub(" ", "_", now() )
time_total <- as.duration(interval(time0, time1))
print(paste0("Total time taken to complete ", numIter, " iterations: ", time_total ))

print("Saving results... ")
write_tsv(results, paste0(data_path, "/IGF1_MR-simulations_iters", numIter,"_", time1,".tsv")) 



