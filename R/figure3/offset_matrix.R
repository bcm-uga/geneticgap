library(LEA)
library(qvalue)
library(gradientForest)

# script to generate data for Figure 3

setwd('~/Documents/Rwork/diversityoffset/')
set.seed(1)

#True positive SNPs
percentage_SNPs_found <- c()
percentage_tp <- c()


# Matrix for each offset
gg_matrix <- c()
rda_matrix <- c()
rda_conf_matrix <- c()
rona_matrix <- c()
rona_conf_matrix <- c()
gf_matrix <- c()
gf_conf_matrix <- c()
logfitnesspred_matrix <- c()
causal_dist_matrix <- c()


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
  # We need genome freq
  Y_freq <- freq_by_pop(Y, pop)
  # We divide Y by 2 in order to be consistent accross all offset
  Y <- Y/2
  X <- cbind(data$var1_step1, data$var2_step1)
  X.pred <- cbind(data$var1_pred, data$var2_pred)
  lfmm.obj <- lfmm2(Y, X, 10, effect.sizes = T)
  pv.obj <- lfmm2.test(lfmm.obj, Y, X, full=T)
  qv.obj <- qvalue(pv.obj$pvalues, fdr.level=0.1)
  snps.set <- which(qv.obj$significant)

  percentage_SNPs_found <- c(percentage_SNPs_found, (sum(data$m2_step1 %in% snps.set) + sum(data$m3_step1 %in% snps.set))/(length(data$m2_step1) + length(data$m3_step1)))
  percentage_tp <- c(percentage_tp, (sum(data$m2_step1 %in% snps.set) + sum(data$m3_step1 %in% snps.set))/length(snps.set))
  
  
  
  # We need Average X variables by pop
  X.pop <- cbind(get_mean_by_pop(data$var1_step1, pop), get_mean_by_pop(data$var2_step1, pop))
  X.pop.pred <- cbind(get_mean_by_pop(data$var1_pred, pop), get_mean_by_pop(data$var2_pred, pop))
  lfmm.obj.freq <- lfmm2(Y_freq, X.pop, 10, effect.sizes = T)
  
  
  
  
  climate_distance <- causal_distance(X.pop, X.pop.pred)
  #RONA
  print("1")
  go.rona <- get_mean_by_pop(go_rona(Y, X, X.pred, snps.set), pop)
  go.rona.conf <- get_mean_by_pop(go_rona_conf(X, X.pred, lfmm.obj, snps.set), pop)
  print("2")
  
  #GENETIC GAP
  gg.obj <- genetic.gap(Y, X, pred.env=X.pred, K=10, scale=T, candidate.loci=snps.set)
  go.gg <- get_mean_by_pop(gg.obj$offset, pop)
  #FITNESS PREd
  fitness_pred <- get_mean_by_pop(data$fitness_pred, pop)
  logfitness <- -log(fitness_pred)
  print("3")
  
  #GF
  gf.obj <- go_gf(Y_freq, X.pop, X.pop.pred, snps.set)
  go.gf <- gf.obj$go
  gf.conf.obj <- go_gf(Y_freq, X.pop, X.pop.pred, snps.set, lfmm.obj.freq@U)
  go.gf.conf <- gf.conf.obj$go
  #RDA
  print("4")
  
  rda.obj <- go_rda(Y, X, X.pred, snps.set)
  go.rda <- get_mean_by_pop(rda.obj, pop)
  rda.obj.conf <- go_rda_conf(Y, X, X.pred,lfmm.obj, snps.set)
  go.rda.conf <- get_mean_by_pop(rda.obj.conf, pop)
  
  
  # Store matrix
  gg_matrix <- rbind(gg_matrix, go.gg)
  
  rda_matrix <- rbind(rda_matrix, go.rda)
  rda_conf_matrix <-  rbind(rda_conf_matrix, go.rda.conf)
  
  rona_matrix <- rbind(rona_matrix, go.rona)
  rona_conf_matrix <- rbind(rona_conf_matrix, go.rona.conf)
  
  gf_matrix <- rbind(gf_matrix, go.gf)
  gf_conf_matrix <- rbind(gf_conf_matrix, go.gf.conf)
  
  logfitnesspred_matrix <- rbind(logfitnesspred_matrix, logfitness)
  
  causal_dist_matrix <- rbind(causal_dist_matrix, climate_distance)
  

}



write.table(gg_matrix, './results/offset_vs_fitness_variation/hchp/gg_matrix')

write.table(rda_matrix, './results/offset_vs_fitness_variation/hchp/rda_matrix')
write.table(rda_conf_matrix, './results/offset_vs_fitness_variation/hchp/rda_conf_matrix')

write.table(rona_matrix, './results/offset_vs_fitness_variation/hchp/rona_matrix')
write.table(rona_conf_matrix, './results/offset_vs_fitness_variation/hchp/rona_conf_matrix')

write.table(gf_matrix, './results/offset_vs_fitness_variation/hchp/gf_matrix')
write.table(gf_conf_matrix, './results/offset_vs_fitness_variation/hchp/gf_conf_matrix')

write.table(logfitnesspred_matrix, './results/offset_vs_fitness_variation/hchp/logfitness_matrix')
write.table(causal_dist_matrix, './results/offset_vs_fitness_variation/hchp/causal_dist_matrix')

write.table(percentage_SNPs_found, './results/offset_vs_fitness_variation/hchp/percentage_snps_found')
write.table(percentage_tp, './results/offset_vs_fitness_variation/hchp/percentage_tp')

#############
# HCWP
#############

#True positive SNPs
percentage_SNPs_found <- c()
percentage_tp <- c()


# Matrix for each offset
gg_matrix <- c()
rda_matrix <- c()
rda_conf_matrix <- c()
rona_matrix <- c()
rona_conf_matrix <- c()
gf_matrix <- c()
gf_conf_matrix <- c()
logfitnesspred_matrix <- c()
causal_dist_matrix <- c()


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
  # We need genome freq
  Y_freq <- freq_by_pop(Y, pop)
  # We divide Y by 2 in order to be consistent accross all offset
  Y <- Y/2
  X <- cbind(data$var1_step1, data$var2_step1)
  X.pred <- cbind(data$var1_pred, data$var2_pred)
  lfmm.obj <- lfmm2(Y, X, 10, effect.sizes = T)
  pv.obj <- lfmm2.test(lfmm.obj, Y, X, full=T)
  qv.obj <- qvalue(pv.obj$pvalues, fdr.level=0.1)
  snps.set <- which(qv.obj$significant)
  
  percentage_SNPs_found <- c(percentage_SNPs_found, (sum(data$m2_step1 %in% snps.set) + sum(data$m3_step1 %in% snps.set))/(length(data$m2_step1) + length(data$m3_step1)))
  percentage_tp <- c(percentage_tp, (sum(data$m2_step1 %in% snps.set) + sum(data$m3_step1 %in% snps.set))/length(snps.set))
  

  # We need Average X variables by pop
  X.pop <- cbind(get_mean_by_pop(data$var1_step1, pop), get_mean_by_pop(data$var2_step1, pop))
  X.pop.pred <- cbind(get_mean_by_pop(data$var1_pred, pop), get_mean_by_pop(data$var2_pred, pop))
  lfmm.obj.freq <- lfmm2(Y_freq, X.pop, 10, effect.sizes = T)
  
  

  
  climate_distance <- causal_distance(X.pop, X.pop.pred)
  #RONA
  go.rona <- get_mean_by_pop(go_rona(Y, X, X.pred, snps.set), pop)
  go.rona.conf <- get_mean_by_pop(go_rona_conf(X, X.pred, lfmm.obj, snps.set), pop)
  #GENETIC GAP
  gg.obj <- genetic.gap(Y, X, pred.env=X.pred, K=10, scale=T, candidate.loci=snps.set)
  go.gg <- get_mean_by_pop(gg.obj$offset, pop)
  #FITNESS PREd
  fitness_pred <- get_mean_by_pop(data$fitness_pred, pop)
  logfitness <- -log(fitness_pred)
  #GF
  gf.obj <- go_gf(Y_freq, X.pop, X.pop.pred, snps.set)
  go.gf <- gf.obj$go
  gf.conf.obj <- go_gf(Y_freq, X.pop, X.pop.pred, snps.set, lfmm.obj.freq@U)
  go.gf.conf <- gf.conf.obj$go
  #RDA
  rda.obj <- go_rda(Y, X, X.pred, snps.set)
  go.rda <- get_mean_by_pop(rda.obj, pop)
  rda.obj.conf <- go_rda_conf(Y, X, X.pred,lfmm.obj, snps.set)
  go.rda.conf <- get_mean_by_pop(rda.obj.conf, pop)
  
  
  # Store matrix
  gg_matrix <- rbind(gg_matrix, go.gg)
  
  rda_matrix <- rbind(rda_matrix, go.rda)
  rda_conf_matrix <-  rbind(rda_conf_matrix, go.rda.conf)
  
  rona_matrix <- rbind(rona_matrix, go.rona)
  rona_conf_matrix <- rbind(rona_conf_matrix, go.rona.conf)
  
  gf_matrix <- rbind(gf_matrix, go.gf)
  gf_conf_matrix <- rbind(gf_conf_matrix, go.gf.conf)
  
  logfitnesspred_matrix <- rbind(logfitnesspred_matrix, logfitness)
  
  causal_dist_matrix <- rbind(causal_dist_matrix, climate_distance)
  
  
}



write.table(gg_matrix, './results/offset_vs_fitness_variation/hcwp/gg_matrix')

write.table(rda_matrix, './results/offset_vs_fitness_variation/hcwp/rda_matrix')
write.table(rda_conf_matrix, './results/offset_vs_fitness_variation/hcwp/rda_conf_matrix')

write.table(rona_matrix, './results/offset_vs_fitness_variation/hcwp/rona_matrix')
write.table(rona_conf_matrix, './results/offset_vs_fitness_variation/hcwp/rona_conf_matrix')

write.table(gf_matrix, './results/offset_vs_fitness_variation/hcwp/gf_matrix')
write.table(gf_conf_matrix, './results/offset_vs_fitness_variation/hcwp/gf_conf_matrix')

write.table(logfitnesspred_matrix, './results/offset_vs_fitness_variation/hcwp/logfitness_matrix')
write.table(causal_dist_matrix, './results/offset_vs_fitness_variation/hcwp/causal_dist_matrix')

write.table(percentage_SNPs_found, './results/offset_vs_fitness_variation/hcwp/percentage_snps_found')
write.table(percentage_tp, './results/offset_vs_fitness_variation/hcwp/percentage_tp')



#############
# LCHP
############

#True positive SNPs
percentage_SNPs_found <- c()
percentage_tp <- c()


# Matrix for each offset
gg_matrix <- c()
rda_matrix <- c()
rda_conf_matrix <- c()
rona_matrix <- c()
rona_conf_matrix <- c()
gf_matrix <- c()
gf_conf_matrix <- c()
logfitnesspred_matrix <- c()
causal_dist_matrix <- c()


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
  # We need genome freq
  Y_freq <- freq_by_pop(Y, pop)
  # We divide Y by 2 in order to be consistent accross all offset
  Y <- Y/2
  X <- cbind(data$var1_step1, data$var2_step1)
  X.pred <- cbind(data$var1_pred, data$var2_pred)
  lfmm.obj <- lfmm2(Y, X, 10, effect.sizes = T)
  pv.obj <- lfmm2.test(lfmm.obj, Y, X, full=T)
  qv.obj <- qvalue(pv.obj$pvalues, fdr.level=0.1)
  snps.set <- which(qv.obj$significant)
  
  percentage_SNPs_found <- c(percentage_SNPs_found, (sum(data$m2_step1 %in% snps.set) + sum(data$m3_step1 %in% snps.set))/(length(data$m2_step1) + length(data$m3_step1)))
  percentage_tp <- c(percentage_tp, (sum(data$m2_step1 %in% snps.set) + sum(data$m3_step1 %in% snps.set))/length(snps.set))
  
  

  # We need Average X variables by pop
  X.pop <- cbind(get_mean_by_pop(data$var1_step1, pop), get_mean_by_pop(data$var2_step1, pop))
  X.pop.pred <- cbind(get_mean_by_pop(data$var1_pred, pop), get_mean_by_pop(data$var2_pred, pop))
  lfmm.obj.freq <- lfmm2(Y_freq, X.pop, 10, effect.sizes = T)
  
  

  
  climate_distance <- causal_distance(X.pop, X.pop.pred)
  #RONA
  go.rona <- get_mean_by_pop(go_rona(Y, X, X.pred, snps.set), pop)
  go.rona.conf <- get_mean_by_pop(go_rona_conf(X, X.pred, lfmm.obj, snps.set), pop)
  #GENETIC GAP
  gg.obj <- genetic.gap(Y, X, pred.env=X.pred, K=10, scale=T, candidate.loci=snps.set)
  go.gg <- get_mean_by_pop(gg.obj$offset, pop)
  #FITNESS PREd
  fitness_pred <- get_mean_by_pop(data$fitness_pred, pop)
  logfitness <- -log(fitness_pred)
  #GF
  gf.obj <- go_gf(Y_freq, X.pop, X.pop.pred, snps.set)
  go.gf <- gf.obj$go
  gf.conf.obj <- go_gf(Y_freq, X.pop, X.pop.pred, snps.set, lfmm.obj.freq@U)
  go.gf.conf <- gf.conf.obj$go
  #RDA
  rda.obj <- go_rda(Y, X, X.pred, snps.set)
  go.rda <- get_mean_by_pop(rda.obj, pop)
  rda.obj.conf <- go_rda_conf(Y, X, X.pred,lfmm.obj, snps.set)
  go.rda.conf <- get_mean_by_pop(rda.obj.conf, pop)
  
  
  # Store matrix
  gg_matrix <- rbind(gg_matrix, go.gg)
  
  rda_matrix <- rbind(rda_matrix, go.rda)
  rda_conf_matrix <-  rbind(rda_conf_matrix, go.rda.conf)
  
  rona_matrix <- rbind(rona_matrix, go.rona)
  rona_conf_matrix <- rbind(rona_conf_matrix, go.rona.conf)
  
  gf_matrix <- rbind(gf_matrix, go.gf)
  gf_conf_matrix <- rbind(gf_conf_matrix, go.gf.conf)
  
  logfitnesspred_matrix <- rbind(logfitnesspred_matrix, logfitness)
  
  causal_dist_matrix <- rbind(causal_dist_matrix, climate_distance)
  
  
}



write.table(gg_matrix, './results/offset_vs_fitness_variation/lchp/gg_matrix')

write.table(rda_matrix, './results/offset_vs_fitness_variation/lchp/rda_matrix')
write.table(rda_conf_matrix, './results/offset_vs_fitness_variation/lchp/rda_conf_matrix')

write.table(rona_matrix, './results/offset_vs_fitness_variation/lchp/rona_matrix')
write.table(rona_conf_matrix, './results/offset_vs_fitness_variation/lchp/rona_conf_matrix')

write.table(gf_matrix, './results/offset_vs_fitness_variation/lchp/gf_matrix')
write.table(gf_conf_matrix, './results/offset_vs_fitness_variation/lchp/gf_conf_matrix')

write.table(logfitnesspred_matrix, './results/offset_vs_fitness_variation/lchp/logfitness_matrix')
write.table(causal_dist_matrix, './results/offset_vs_fitness_variation/lchp/causal_dist_matrix')

write.table(percentage_SNPs_found, './results/offset_vs_fitness_variation/lchp/percentage_snps_found')
write.table(percentage_tp, './results/offset_vs_fitness_variation/lchp/percentage_tp')


###################
#   LCWP
###################

#True positive SNPs
percentage_SNPs_found <- c()
percentage_tp <- c()


# Matrix for each offset
gg_matrix <- c()
rda_matrix <- c()
rda_conf_matrix <- c()
rona_matrix <- c()
rona_conf_matrix <- c()
gf_matrix <- c()
gf_conf_matrix <- c()
logfitnesspred_matrix <- c()
causal_dist_matrix <- c()


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
  # We need genome freq
  Y_freq <- freq_by_pop(Y, pop)
  # We divide Y by 2 in order to be consistent accross all offset
  Y <- Y/2
  X <- cbind(data$var1_step1, data$var2_step1)
  X.pred <- cbind(data$var1_pred, data$var2_pred)
  lfmm.obj <- lfmm2(Y, X, 10, effect.sizes = T)
  pv.obj <- lfmm2.test(lfmm.obj, Y, X, full=T)
  qv.obj <- qvalue(pv.obj$pvalues, fdr.level=0.1)
  snps.set <- which(qv.obj$significant)
  
  percentage_SNPs_found <- c(percentage_SNPs_found, (sum(data$m2_step1 %in% snps.set) + sum(data$m3_step1 %in% snps.set))/(length(data$m2_step1) + length(data$m3_step1)))
  percentage_tp <- c(percentage_tp, (sum(data$m2_step1 %in% snps.set) + sum(data$m3_step1 %in% snps.set))/length(snps.set))
  
  

  # We need Average X variables by pop
  X.pop <- cbind(get_mean_by_pop(data$var1_step1, pop), get_mean_by_pop(data$var2_step1, pop))
  X.pop.pred <- cbind(get_mean_by_pop(data$var1_pred, pop), get_mean_by_pop(data$var2_pred, pop))
  lfmm.obj.freq <- lfmm2(Y_freq, X.pop, 10, effect.sizes = T)
  
  

  
  climate_distance <- causal_distance(X.pop, X.pop.pred)
  #RONA
  go.rona <- get_mean_by_pop(go_rona(Y, X, X.pred, snps.set), pop)
  go.rona.conf <- get_mean_by_pop(go_rona_conf(X, X.pred, lfmm.obj, snps.set), pop)
  #GENETIC GAP
  gg.obj <- genetic.gap(Y, X, pred.env=X.pred, K=10, scale=T, candidate.loci=snps.set)
  go.gg <- get_mean_by_pop(gg.obj$offset, pop)
  #FITNESS PREd
  fitness_pred <- get_mean_by_pop(data$fitness_pred, pop)
  logfitness <- -log(fitness_pred)
  #GF
  gf.obj <- go_gf(Y_freq, X.pop, X.pop.pred, snps.set)
  go.gf <- gf.obj$go
  gf.conf.obj <- go_gf(Y_freq, X.pop, X.pop.pred, snps.set, lfmm.obj.freq@U)
  go.gf.conf <- gf.conf.obj$go
  #RDA
  rda.obj <- go_rda(Y, X, X.pred, snps.set)
  go.rda <- get_mean_by_pop(rda.obj, pop)
  rda.obj.conf <- go_rda_conf(Y, X, X.pred,lfmm.obj, snps.set)
  go.rda.conf <- get_mean_by_pop(rda.obj.conf, pop)
  
  
  # Store matrix
  gg_matrix <- rbind(gg_matrix, go.gg)
  
  rda_matrix <- rbind(rda_matrix, go.rda)
  rda_conf_matrix <-  rbind(rda_conf_matrix, go.rda.conf)
  
  rona_matrix <- rbind(rona_matrix, go.rona)
  rona_conf_matrix <- rbind(rona_conf_matrix, go.rona.conf)
  
  gf_matrix <- rbind(gf_matrix, go.gf)
  gf_conf_matrix <- rbind(gf_conf_matrix, go.gf.conf)
  
  logfitnesspred_matrix <- rbind(logfitnesspred_matrix, logfitness)
  
  causal_dist_matrix <- rbind(causal_dist_matrix, climate_distance)
  
  
}



write.table(gg_matrix, './results/offset_vs_fitness_variation/lcwp/gg_matrix')

write.table(rda_matrix, './results/offset_vs_fitness_variation/lcwp/rda_matrix')
write.table(rda_conf_matrix, './results/offset_vs_fitness_variation/lcwp/rda_conf_matrix')

write.table(rona_matrix, './results/offset_vs_fitness_variation/lcwp/rona_matrix')
write.table(rona_conf_matrix, './results/offset_vs_fitness_variation/lcwp/rona_conf_matrix')

write.table(gf_matrix, './results/offset_vs_fitness_variation/lcwp/gf_matrix')
write.table(gf_conf_matrix, './results/offset_vs_fitness_variation/lcwp/gf_conf_matrix')

write.table(logfitnesspred_matrix, './results/offset_vs_fitness_variation/lcwp/logfitness_matrix')
write.table(causal_dist_matrix, './results/offset_vs_fitness_variation/lcwp/causal_dist_matrix')

write.table(percentage_SNPs_found, './results/offset_vs_fitness_variation/lcwp/percentage_snps_found')
write.table(percentage_tp, './results/offset_vs_fitness_variation/lcwp/percentage_tp')


