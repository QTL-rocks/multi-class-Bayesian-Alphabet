setwd("/Users/ziguiwang/Documents/Animal/Functional_Annoation/FAANG_simulation/reports/MSU/")

marker_var_no_group = read.csv("no_snp_group/BayesC_result/MCMC_samples_marker_effects_variances_genotype1.txt",
                               header = FALSE)

mean(marker_var_no_group$V1)


marker_var_group1 = read.csv("no_snp_group/BayesC_result1/MCMC_samples_marker_effects_variances_genotype1.txt",
                             header = FALSE)

group_number = 1

path = paste("snp_group/BayesC_result1/MCMC_samples_marker_effects_variances_genotype", 
      as.character(group_number),".txt", sep = "")


marker_var_group1 = read.csv(path,header = FALSE)

read_marker_effect_var <-function(group_number){
  path = paste("snp_group/BayesC_result/MCMC_samples_marker_effects_variances_genotype", 
               as.character(group_number),".txt", sep = "")
  
  
  marker_var_group1 = read.csv(path,header = FALSE)
  
  return(mean(marker_var_group1$V1 ) )
}
read_marker_effect_var(1)
for (i in 1:5){
  print(read_marker_effect_var(i))
}

