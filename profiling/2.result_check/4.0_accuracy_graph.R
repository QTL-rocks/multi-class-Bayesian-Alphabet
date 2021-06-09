### This R script is used to generate the GWAS graph for single trait analysis in paper "Single-trait and 
###  Multiple-trait Genomic Prediction from Multi-class Bayesian Alphabet Models using Biological Information"

## Since the graph generation process for single-trait result and multi-trait result are same, here we only show 
## how to generate the result graph for single-trait result.

# load packages
library(greekLetters)  
library(ggplot2)
library(cowplot)

# set the working directory
setwd("/Users/ziguiwang/Documents/Animal/Functional_Annoation/FAANG_simulation/reports/accuracy/MSU/")

################################################################################
# Function setup
################################################################################

# This function is used to generate the accuracy scatter plot
plot_function<-function(df,significant_ans){
  lower_bound = 0.3
  upper_bound = 0.7
    result_plot = ggplot(df, aes(x=standard, y= Multi,color = pheno)) + 
      geom_point(alpha=0.8) +
      xlim(lower_bound, upper_bound) +
      ylim(lower_bound, upper_bound) +
      xlab("conventional")+
      ylab("multi-class") +
      geom_abline(intercept = 0 , slope = 1) +
      theme(legend.position = "none") 
  
  return(result_plot)
}

# this function is used to transform the result dataframe to the dataframe that can be used in ggplot2
get_method_df<-function(method_name_x,method_name_y,standard_all_accuracy,Multi_all_accuracy){
  
  # find the standard and MUlti accuracy for this method
  method_standard_accuracy = standard_all_accuracy$accuracy[standard_all_accuracy$method == method_name_x]
  method_Multi_accuracy = Multi_all_accuracy$accuracy[Multi_all_accuracy$method == method_name_y]
  
  rep = standard_all_accuracy$rep[standard_all_accuracy$method == method_name_x]
  pheno = standard_all_accuracy$pheno[standard_all_accuracy$method == method_name_x]
  
  # define the data frame for ggplot
  method_df = data.frame("standard" = method_standard_accuracy ,
                         "Multi" = method_Multi_accuracy,method = method_name_x, rep = rep,pheno = as.character(pheno))
  
  return(method_df)
}

################################################################################
# Read the result
################################################################################

# set the initial result data frame
standard_all_accuracy = data.frame("method" = c(),"accuracy" = c(),rep = c())
Multi_all_accuracy = data.frame("method" = c(),"accuracy" = c(),rep = c())


# loop through all the results
for (i in 1:5){          # loop through 5 replications
for (j in 1:30){         # loop through 30 simulated datasets
  
  # read the conventional Bayesian Alphabet result
  standard_accuracy = read.csv(paste("no_snp_group/rep_",as.character(i) ,"_phenotype_",as.character(j),".csv",sep = ""))
  # read the multi-class Bayesian Alphabet result
  Multi_accuracy = read.csv(paste("snp_group/rep_", as.character(i),"_phenotype_",as.character(j),".csv",sep = ""))
  
  # add the replication and phenotype index to the data frame
  standard_accuracy$rep = i
  Multi_accuracy$rep = i
  standard_accuracy$phenotype = j
  Multi_accuracy$phenotype = j
  
  # merge the current result to the final result
  standard_all_accuracy = rbind(standard_all_accuracy,standard_accuracy)
  Multi_all_accuracy = rbind(Multi_all_accuracy,Multi_accuracy)
}
}




################################################################################
# genreate the graph
################################################################################

# Get the dataframe for ggplot2
RRBLUP_df = get_method_df("RR-BLUP","RR-BLUP",standard_all_accuracy,Multi_all_accuracy)
BayesA_df = get_method_df("BayesA","BayesA",standard_all_accuracy,Multi_all_accuracy)
BayesB_df = get_method_df("BayesB","BayesB",standard_all_accuracy,Multi_all_accuracy)
BayesC_df = get_method_df("BayesC","BayesC",standard_all_accuracy,Multi_all_accuracy)
BayesL_df = get_method_df("BayesL","BayesL",standard_all_accuracy,Multi_all_accuracy)
ensemble_df = get_method_df("ensemble","ensemble",standard_all_accuracy,Multi_all_accuracy)


# Do the paired t test with 0.1 significance level
t.test(BayesA_df$standard ,BayesA_df$Multi,conf.level = 0.9) 
t.test(BayesB_df$standard ,BayesB_df$Multi,conf.level = 0.9) 
t.test(BayesC_df$standard ,BayesC_df$Multi,conf.level = 0.9) 
t.test(BayesL_df$standard ,BayesL_df$Multi,conf.level = 0.9) 
t.test(RRBLUP_df$standard ,RRBLUP_df$Multi,conf.level = 0.9) 
t.test(ensemble_df$standard ,ensemble_df$Multi,conf.level = 0.9) 


# merge the results from 6 methods
overall_df = rbind(RRBLUP_df,BayesA_df,BayesB_df,BayesC_df,BayesL_df,ensemble_df)
overall_df[overall_df$method == "BayesL",3] = "Bayesian LASSO"
overall_df[overall_df$method == "BayesC",3] = paste("BayesC",greeks("pi"),sep = "")


overall_df$method_f = factor(overall_df$method, levels=c('RR-BLUP','BayesA',
                                                         'BayesB',paste("BayesC",greeks("pi"),sep = ""),
                                                         "Bayesian LASSO","ensemble"))

# Generate the scatter plot
ggplot(overall_df , aes(x=standard, y= Multi, color = pheno)) +
  geom_abline(intercept = 0 , slope = 1)+
  geom_point() +
  facet_wrap(~ method_f) +
  guides(fill = FALSE) +  # to remove the legend
  theme_bw() + 
  guides(fill = FALSE) + 
  xlab("conventional Bayesian Alphabet")+
  ylab("multi-class Bayesian Alphabet")+
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        #axis.text.x=element_blank(),
        #axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        #axis.text.y=element_blank(),
        #axis.ticks.y=element_blank(),
        strip.text = element_text(size = 20)) 


