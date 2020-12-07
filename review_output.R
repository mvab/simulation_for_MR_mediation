library(dplyr)
library(readr)
library(tibble)


file <- "data/age_at_menopause_MR-simulations_iters1000_2020-12-07_00-50-56.tsv"
data <- read_tsv(file) 
colnames(data)

means <- data %>% select(-starts_with("param")) %>% 
   summarise(across(everything(), mean))

sds <- data %>% select(-starts_with("param")) %>% 
  summarise(across(everything(), sd))

out <- rbind(means, sds) %>% 
  t() %>% as.data.frame() %>% 
  rownames_to_column('measure')
colnames(out)[2:3] <- c("mean", "sd")

out %>%  write_tsv(paste0(strsplit(file, "MR")[[1]][1], "summary.tsv"))
out %>% print()

#     measure                          mean           sd
# EO_total_beta                -0.395213986 0.0111551336
# EM_total_beta                -0.784233588 0.0127446164
# MO_total_beta                -0.018770618 0.0122625478
# EO_direct_beta               -0.400977382 0.0147544055
# MO_direct_beta               -0.011993242 0.0130537477
# EO_total_se                   0.010825020 0.0007964317
# EM_total_se                   0.012810320 0.0009622165
# MO_total_se                   0.011985474 0.0010729079
# EO_direct_se                  0.014927585 0.0009012914
# MO_direct_se                  0.013184216 0.0008237578
# indirect_b_difference         0.005763396 0.0102506173
# indirect_se_difference_PoE    0.018444721 0.0011191846
# indirect_b_product_v1         0.014717273 0.0096130378
# indirect_se_product_v1_PoE    0.017570926 0.0010461605
# indirect_se_product_v1_delta  0.009404349 0.0008604131
# indirect_b_product_v2         0.009405936 0.0102391091
# indirect_se_product_v2_PoE    0.018403531 0.0009173742
# indirect_se_product_v2_delta  0.010342437 0.0006746829




# What we are looking for is whether the mean of the SE calculated using either the PoE or Delta methods 
# are similar to the standard deviation of the relevant effect estimate across all of the simulations. 

out %>% 
  filter(grepl( "indirect", measure)) %>% 
  mutate(mean = ifelse(grepl("b_", measure), NA, mean),
         sd = ifelse(grepl("se_", measure), NA, sd))



#                       measure        mean          sd
# 1        indirect_b_difference          NA 0.002883507
# 2   indirect_se_difference_PoE 0.015411137          NA
# 3        indirect_b_product_v1          NA 0.001833279
# 4   indirect_se_product_v1_PoE 0.014903687          NA
# 5 indirect_se_product_v1_delta 0.001828008          NA
# 6        indirect_b_product_v2          NA 0.001902260
# 7   indirect_se_product_v2_PoE 0.015120794          NA
# 8 indirect_se_product_v2_delta 0.001943369          NA
