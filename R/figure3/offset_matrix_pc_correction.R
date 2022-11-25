library(LEA)
library(qvalue)
library(gradientForest)

# script to generate data for Figure 3

setwd('~/Documents/Rwork/diversityoffset/')
set.seed(1)



# Matrix for each offset
rda_conf_matrix <- c()
rona_conf_matrix <- c()
gf_conf_matrix <- c()

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
  Y_freq <- freq_by_pop(Y, pop)
  Y <- Y/2
  pc.obj <- prcomp(Y)
  X <- cbind(data$var1_step1, data$var2_step1)
  X.pred <- cbind(data$var1_pred, data$var2_pred)
  lfmm.obj <- lfmm2(Y, X, 10, effect.sizes = T)
  pv.obj <- lfmm2.test(lfmm.obj, Y, X, full=T)
  qv.obj <- qvalue(pv.obj$pvalues, fdr.level=0.1)
  snps.set <- which(qv.obj$significant)


  # We need Average X variables by pop
  X.pop <- cbind(get_mean_by_pop(data$var1_step1, pop), get_mean_by_pop(data$var2_step1, pop))
  X.pop.pred <- cbind(get_mean_by_pop(data$var1_pred, pop), get_mean_by_pop(data$var2_pred, pop))
  lfmm.obj.freq <- lfmm2(Y_freq, X.pop, 10, effect.sizes = T)
  pc.freq <- prcomp(Y_freq)
  

  #RONA
  go.rona.conf <- get_mean_by_pop(go_rona_pc(as.matrix(Y), as.matrix(X), as.matrix(X.pred), pc.obj$x[,1:10], snps.set), pop)
  
  
  #GF

  gf.conf.obj <- go_gf(Y_freq, X.pop, X.pop.pred, snps.set, pc.freq$x[,1:10])
  go.gf.conf <- gf.conf.obj$go
  #RDA

  rda.obj.conf <- go_rda_pc(as.matrix(Y), as.matrix(X), as.matrix(X.pred), pc.obj$x[,1:10], snps.set)
  go.rda.conf <- get_mean_by_pop(rda.obj.conf, pop)
  
  
  # Store matrix

  rda_conf_matrix <-  rbind(rda_conf_matrix, go.rda.conf)
  rona_conf_matrix <- rbind(rona_conf_matrix, go.rona.conf)
  gf_conf_matrix <- rbind(gf_conf_matrix, go.gf.conf)

}




write.table(rda_conf_matrix, './results/offset_vs_fitness_variation/hchp/rda_pc_matrix')
write.table(rona_conf_matrix, './results/offset_vs_fitness_variation/hchp/rona_pc_matrix')
write.table(gf_conf_matrix, './results/offset_vs_fitness_variation/hchp/gf_pc_matrix')



#############
# HCWP
#############


# Matrix for each offset
rda_conf_matrix <- c()
rona_conf_matrix <- c()
gf_conf_matrix <- c()

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
  Y_freq <- freq_by_pop(Y, pop)
  Y <- Y/2
  pc.obj <- prcomp(Y)
  X <- cbind(data$var1_step1, data$var2_step1)
  X.pred <- cbind(data$var1_pred, data$var2_pred)
  lfmm.obj <- lfmm2(Y, X, 10, effect.sizes = T)
  pv.obj <- lfmm2.test(lfmm.obj, Y, X, full=T)
  qv.obj <- qvalue(pv.obj$pvalues, fdr.level=0.1)
  snps.set <- which(qv.obj$significant)
  
  
  # We need Average X variables by pop
  X.pop <- cbind(get_mean_by_pop(data$var1_step1, pop), get_mean_by_pop(data$var2_step1, pop))
  X.pop.pred <- cbind(get_mean_by_pop(data$var1_pred, pop), get_mean_by_pop(data$var2_pred, pop))
  lfmm.obj.freq <- lfmm2(Y_freq, X.pop, 10, effect.sizes = T)
  pc.freq <- prcomp(Y_freq)
  
  
  #RONA
  go.rona.conf <- get_mean_by_pop(go_rona_pc(as.matrix(Y), as.matrix(X), as.matrix(X.pred), pc.obj$x[,1:10], snps.set), pop)
  
  
  #GF
  
  gf.conf.obj <- go_gf(Y_freq, X.pop, X.pop.pred, snps.set, pc.freq$x[,1:10])
  go.gf.conf <- gf.conf.obj$go
  #RDA
  
  rda.obj.conf <- go_rda_pc(as.matrix(Y),as.matrix(X), as.matrix(X.pred), pc.obj$x[,1:10], snps.set)
  go.rda.conf <- get_mean_by_pop(rda.obj.conf, pop)
  
  
  # Store matrix
  
  rda_conf_matrix <-  rbind(rda_conf_matrix, go.rda.conf)
  rona_conf_matrix <- rbind(rona_conf_matrix, go.rona.conf)
  gf_conf_matrix <- rbind(gf_conf_matrix, go.gf.conf)
  
}




write.table(rda_conf_matrix, './results/offset_vs_fitness_variation/hcwp/rda_pc_matrix')
write.table(rona_conf_matrix, './results/offset_vs_fitness_variation/hcwp/rona_pc_matrix')
write.table(gf_conf_matrix, './results/offset_vs_fitness_variation/hcwp/gf_pc_matrix')


#############
# LCHP
############


# Matrix for each offset
rda_conf_matrix <- c()
rona_conf_matrix <- c()
gf_conf_matrix <- c()

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
  Y_freq <- freq_by_pop(Y, pop)
  Y <- Y/2
  pc.obj <- prcomp(Y)
  X <- cbind(data$var1_step1, data$var2_step1)
  X.pred <- cbind(data$var1_pred, data$var2_pred)
  lfmm.obj <- lfmm2(Y, X, 10, effect.sizes = T)
  pv.obj <- lfmm2.test(lfmm.obj, Y, X, full=T)
  qv.obj <- qvalue(pv.obj$pvalues, fdr.level=0.1)
  snps.set <- which(qv.obj$significant)
  
  
  # We need Average X variables by pop
  X.pop <- cbind(get_mean_by_pop(data$var1_step1, pop), get_mean_by_pop(data$var2_step1, pop))
  X.pop.pred <- cbind(get_mean_by_pop(data$var1_pred, pop), get_mean_by_pop(data$var2_pred, pop))
  lfmm.obj.freq <- lfmm2(Y_freq, X.pop, 10, effect.sizes = T)
  pc.freq <- prcomp(Y_freq)
  
  
  #RONA
  go.rona.conf <- get_mean_by_pop(go_rona_pc(as.matrix(Y), as.matrix(X), as.matrix(X.pred), pc.obj$x[,1:10], snps.set), pop)
  
  
  #GF
  
  gf.conf.obj <- go_gf(Y_freq, X.pop, X.pop.pred, snps.set, pc.freq$x[,1:10])
  go.gf.conf <- gf.conf.obj$go
  #RDA
  
  rda.obj.conf <- go_rda_pc(as.matrix(Y), as.matrix(X), as.matrix(X.pred), pc.obj$x[,1:10], snps.set)
  go.rda.conf <- get_mean_by_pop(rda.obj.conf, pop)
  
  
  # Store matrix
  
  rda_conf_matrix <-  rbind(rda_conf_matrix, go.rda.conf)
  rona_conf_matrix <- rbind(rona_conf_matrix, go.rona.conf)
  gf_conf_matrix <- rbind(gf_conf_matrix, go.gf.conf)
  
}




write.table(rda_conf_matrix, './results/offset_vs_fitness_variation/lchp/rda_pc_matrix')
write.table(rona_conf_matrix, './results/offset_vs_fitness_variation/lchp/rona_pc_matrix')
write.table(gf_conf_matrix, './results/offset_vs_fitness_variation/lchp/gf_pc_matrix')


###################
#   LCWP
###################



# Matrix for each offset
rda_conf_matrix <- c()
rona_conf_matrix <- c()
gf_conf_matrix <- c()

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
  Y_freq <- freq_by_pop(Y, pop)
  Y <- Y/2
  pc.obj <- prcomp(Y)
  X <- cbind(data$var1_step1, data$var2_step1)
  X.pred <- cbind(data$var1_pred, data$var2_pred)
  lfmm.obj <- lfmm2(Y, X, 10, effect.sizes = T)
  pv.obj <- lfmm2.test(lfmm.obj, Y, X, full=T)
  qv.obj <- qvalue(pv.obj$pvalues, fdr.level=0.1)
  snps.set <- which(qv.obj$significant)
  
  
  # We need Average X variables by pop
  X.pop <- cbind(get_mean_by_pop(data$var1_step1, pop), get_mean_by_pop(data$var2_step1, pop))
  X.pop.pred <- cbind(get_mean_by_pop(data$var1_pred, pop), get_mean_by_pop(data$var2_pred, pop))
  lfmm.obj.freq <- lfmm2(Y_freq, X.pop, 10, effect.sizes = T)
  pc.freq <- prcomp(Y_freq)
  
  
  #RONA
  go.rona.conf <- get_mean_by_pop(go_rona_pc(as.matrix(Y), as.matrix(X), as.matrix(X.pred), pc.obj$x[,1:10], snps.set), pop)
  
  
  #GF
  
  gf.conf.obj <- go_gf(Y_freq, X.pop, X.pop.pred, snps.set, pc.freq$x[,1:10])
  go.gf.conf <- gf.conf.obj$go
  #RDA
  
  rda.obj.conf <- go_rda_pc(as.matrix(Y), as.matrix(X), as.matrix(X.pred), pc.obj$x[,1:10], snps.set)
  go.rda.conf <- get_mean_by_pop(rda.obj.conf, pop)
  
  
  # Store matrix
  
  rda_conf_matrix <-  rbind(rda_conf_matrix, go.rda.conf)
  rona_conf_matrix <- rbind(rona_conf_matrix, go.rona.conf)
  gf_conf_matrix <- rbind(gf_conf_matrix, go.gf.conf)
  
}




write.table(rda_conf_matrix, './results/offset_vs_fitness_variation/lcwp/rda_pc_matrix')
write.table(rona_conf_matrix, './results/offset_vs_fitness_variation/lcwp/rona_pc_matrix')
write.table(gf_conf_matrix, './results/offset_vs_fitness_variation/lcwp/gf_pc_matrix')

