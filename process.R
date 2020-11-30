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



simulate_results <- function(n_iter, snps_m, beta_em, beta_mo){
  
  #data generation
  MR.dat <- datagen(snps_m, beta_em, beta_mo)

  #analysis
  ex.dat <- MR.dat[1:100,]
  med.dat <- MR.dat[101:100+snps_m,]
  
  #IVW estimation - exposure -> outcome 
  eo <- lm(ex.dat$Out_beta~-1+ex.dat$Ex_beta, weights = 1/ex.dat$Out_se^2)
  eo_value <- eo$coefficients %>% as.numeric()
  
  #IVW estimation - exposure -> mediator
  em <- lm(ex.dat$Med_beta~-1+ex.dat$Ex_beta, weights = 1/ex.dat$Med_se^2)
  em_value <- em$coefficients %>% as.numeric()
  
  #IVW estimation - mediator -> outcome
  mo <- lm(med.dat$Out_beta~-1+med.dat$Med_beta, weights = 1/med.dat$Out_se^2)
  mo_value <- mo$coefficients %>% as.numeric()
  
  #IVW - MVMR estimation
  emo <- lm(MR.dat$Out_beta~-1+MR.dat$Ex_beta+MR.dat$Med_beta, weights = 1/MR.dat$Out_se^2)
  emo_value <- emo$coefficients %>% as.numeric()
  
  output_vector <-  data.frame(eo_value, em_value, mo_value, emo_value[1], emo_value[2])  

  output_df<-rbind(data.frame(), output_vector)
  
  print(paste0("Simulation iteration number ", n_iter, "completed at:", gsub(" ", "_", Sys.time())  ))
  
  return(output_df)
}


n_iter = seq(1:numIter)
results_list <- mclapply(n_iter, simulate_results, mc.cores = numCores,
                                snps_m = snps_m,
                                beta_em = beta_em, 
                                beta_mo = beta_mo)


# convert list of vectors to df
results <- bind_rows(lapply(results_list, as.data.frame.list))
colnames(results) <- c("EO_total", "EM_total", "MO_total", "EO_direct", 'MO_direct')

print("Saving results... ")
time1 <- gsub(" ", "_", now())
write_tsv(results, paste0(data_path, "/IGF1_MR-simulations_iters", numIter,"_", time1,".tsv")) 


time_total <- interval(time0, time1)
print(paste0("Total time taken to complete ", n_iter, " iterations: ", time_total ))


