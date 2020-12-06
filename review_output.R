library(dplyr)
library(readr)


data <- read_tsv("data/IGF1_MR-simulations_iters1000_2020-12-04_12:09:22.tsv") 
colnames(data)

means <- data %>% select(-starts_with("param")) %>% 
   summarise(across(everything(), mean))

sds <- data %>% select(-starts_with("param")) %>% 
  summarise(across(everything(), sd))

out <- rbind(means, sds) %>% t()
colnames(out) <- c("mean", "sd")
