library(LEA)
library(qvalue)
library(gradientForest)

# script to generate data for Figure 3

setwd('~/Documents/Rwork/diversityoffset/')
set.seed(1)



# Matrix for each offset
causaldist_matrix <- c()


# We iterate over all simulation scenario 
for (i in seq(1,30)){
  
  print("=======")
  print(i)
  print("=======")
  directory_name <-  paste('./data/expfit_2variables_poly/', as.character(i), sep="")
  directory_name <-  paste(directory_name, "/", sep="")
  data <- extract_data(directory_name, MAF=0.01)
  pop <- get_pop(data$pos_step1, 12, 12, 12, 12)
  L <- dim(data$genome_step1)[2]
  
  # Get snps.set
  
  Y <- data$genome_step1
  X <- cbind(data$var1_step1, data$var2_step1)
  X.pred <- cbind(data$var1_pred, data$var2_pred)
  
  # We need Average X variables by pop
  X.pop <- cbind(get_mean_by_pop(data$var1_step1, pop), get_mean_by_pop(data$var2_step1, pop))
  X.pop.pred <- cbind(get_mean_by_pop(data$var1_pred, pop), get_mean_by_pop(data$var2_pred, pop))
  
  climate_distance <- causal_distance(X.pop, X.pop.pred)

  
  # Store matrix
  causaldist_matrix <- rbind(causaldist_matrix, climate_distance)
  
}



write.table(causaldist_matrix, './results/offset_vs_fitness_variation/hchp/causaldist_matrix')


# Matrix for each offset
causaldist_matrix <- c()


# We iterate over all simulation scenario 
for (i in seq(1,30)){
  
  print("=======")
  print(i)
  print("=======")
  directory_name <-  paste('./data/expfit_2variables/', as.character(i), sep="")
  directory_name <-  paste(directory_name, "/", sep="")
  data <- extract_data(directory_name, MAF=0.01)
  pop <- get_pop(data$pos_step1, 12, 12, 12, 12)
  L <- dim(data$genome_step1)[2]
  
  # Get snps.set
  
  Y <- data$genome_step1
  X <- cbind(data$var1_step1, data$var2_step1)
  X.pred <- cbind(data$var1_pred, data$var2_pred)
  
  # We need Average X variables by pop
  X.pop <- cbind(get_mean_by_pop(data$var1_step1, pop), get_mean_by_pop(data$var2_step1, pop))
  X.pop.pred <- cbind(get_mean_by_pop(data$var1_pred, pop), get_mean_by_pop(data$var2_pred, pop))
  
  climate_distance <- causal_distance(X.pop, X.pop.pred)
  
  
  # Store matrix
  causaldist_matrix <- rbind(causaldist_matrix, climate_distance)
  
}



write.table(causaldist_matrix, './results/offset_vs_fitness_variation/hcwp/causaldist_matrix')


# Matrix for each offset
causaldist_matrix <- c()


# We iterate over all simulation scenario 
for (i in seq(1,30)){
  
  print("=======")
  print(i)
  print("=======")
  directory_name <-  paste('./data/poly_exp/', as.character(i), sep="")
  directory_name <-  paste(directory_name, "/", sep="")
  data <- extract_data(directory_name, MAF=0.01)
  pop <- get_pop(data$pos_step1, 12, 12, 12, 12)
  L <- dim(data$genome_step1)[2]
  
  # Get snps.set
  
  Y <- data$genome_step1
  X <- cbind(data$var1_step1, data$var2_step1)
  X.pred <- cbind(data$var1_pred, data$var2_pred)
  
  # We need Average X variables by pop
  X.pop <- cbind(get_mean_by_pop(data$var1_step1, pop), get_mean_by_pop(data$var2_step1, pop))
  X.pop.pred <- cbind(get_mean_by_pop(data$var1_pred, pop), get_mean_by_pop(data$var2_pred, pop))
  
  climate_distance <- causal_distance(X.pop, X.pop.pred)
  
  
  # Store matrix
  causaldist_matrix <- rbind(causaldist_matrix, climate_distance)
  
}



write.table(causaldist_matrix, './results/offset_vs_fitness_variation/lchp/causaldist_matrix')


# Matrix for each offset
causaldist_matrix <- c()


# We iterate over all simulation scenario 
for (i in seq(1,30)){
  
  print("=======")
  print(i)
  print("=======")
  directory_name <-  paste('./data/poly_small/', as.character(i), sep="")
  directory_name <-  paste(directory_name, "/", sep="")
  data <- extract_data(directory_name, MAF=0.01)
  pop <- get_pop(data$pos_step1, 12, 12, 12, 12)
  L <- dim(data$genome_step1)[2]
  
  # Get snps.set
  
  Y <- data$genome_step1
  X <- cbind(data$var1_step1, data$var2_step1)
  X.pred <- cbind(data$var1_pred, data$var2_pred)
  
  # We need Average X variables by pop
  X.pop <- cbind(get_mean_by_pop(data$var1_step1, pop), get_mean_by_pop(data$var2_step1, pop))
  X.pop.pred <- cbind(get_mean_by_pop(data$var1_pred, pop), get_mean_by_pop(data$var2_pred, pop))
  
  climate_distance <- causal_distance(X.pop, X.pop.pred)
  
  
  # Store matrix
  causaldist_matrix <- rbind(causaldist_matrix, climate_distance)
  
}



write.table(causaldist_matrix, './results/offset_vs_fitness_variation/lcwp/causaldist_matrix')