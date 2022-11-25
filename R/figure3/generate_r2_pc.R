# The aim here is to generate R2 for the different scenario


setwd('~/Documents/Rwork/diversityoffset/')
set.seed(1)


# HCWP

matrix_rda_conf <- read.table('./results/offset_vs_fitness_variation/hcwp/rda_pc_matrix', row.names = NULL)[,2:101]
matrix_rona_conf <- read.table('./results/offset_vs_fitness_variation/hcwp/rona_pc_matrix', row.names = NULL)[,2:101]
matrix_gf_conf <- read.table('./results/offset_vs_fitness_variation/hcwp/gf_pc_matrix', row.names = NULL)[,2:101]
matrix_fitness <- read.table('./results/offset_vs_fitness_variation/hcwp/logfitness_matrix', row.names = NULL)[,2:101]


rda_conf_r2 <- c()
rona_conf_r2 <- c()
gf_conf_r2 <- c()

for (i in seq(1,30)){
  rda_conf_r2 <- c(rda_conf_r2, summary(lm(as.numeric(matrix_fitness[i,]) ~ as.numeric(matrix_rda_conf[i,])^2))$adj.r.squared)
  rona_conf_r2 <- c(rona_conf_r2, summary(lm(as.numeric(matrix_fitness[i,]) ~ as.numeric(matrix_rona_conf[i,])^2))$adj.r.squared)
  gf_conf_r2 <- c(gf_conf_r2, summary(lm(as.numeric(matrix_fitness[i,]) ~ as.numeric(matrix_gf_conf[i,])^2))$adj.r.squared)
}



result_r2_hcwp <- cbind( rda_conf_r2,  rona_conf_r2,  gf_conf_r2)
colnames(result_r2_hcwp) <- c( 'rda_pc', 'rona_pc', 'gf_pc')
write.table(result_r2_hcwp, './results/offset_vs_fitness_variation/hcwp/r2_pc_matrix')



# HCHP
matrix_rda_conf <- read.table('./results/offset_vs_fitness_variation/hchp/rda_pc_matrix', row.names = NULL)[,2:101]
matrix_rona_conf <- read.table('./results/offset_vs_fitness_variation/hchp/rona_pc_matrix', row.names = NULL)[,2:101]
matrix_gf_conf <- read.table('./results/offset_vs_fitness_variation/hchp/gf_pc_matrix', row.names = NULL)[,2:101]
matrix_fitness <- read.table('./results/offset_vs_fitness_variation/hchp/logfitness_matrix', row.names = NULL)[,2:101]


rda_conf_r2 <- c()
rona_conf_r2 <- c()
gf_conf_r2 <- c()

for (i in seq(1,30)){
  rda_conf_r2 <- c(rda_conf_r2, summary(lm(as.numeric(matrix_fitness[i,]) ~ as.numeric(matrix_rda_conf[i,])^2))$adj.r.squared)
  rona_conf_r2 <- c(rona_conf_r2, summary(lm(as.numeric(matrix_fitness[i,]) ~ as.numeric(matrix_rona_conf[i,])^2))$adj.r.squared)
  gf_conf_r2 <- c(gf_conf_r2, summary(lm(as.numeric(matrix_fitness[i,]) ~ as.numeric(matrix_gf_conf[i,])^2))$adj.r.squared)
}



result_r2_hcwp <- cbind( rda_conf_r2,  rona_conf_r2,  gf_conf_r2)
colnames(result_r2_hcwp) <- c( 'rda_pc', 'rona_pc', 'gf_pc')
write.table(result_r2_hcwp, './results/offset_vs_fitness_variation/hchp/r2_pc_matrix')

# LCWP
matrix_rda_conf <- read.table('./results/offset_vs_fitness_variation/lcwp/rda_pc_matrix', row.names = NULL)[,2:101]
matrix_rona_conf <- read.table('./results/offset_vs_fitness_variation/lcwp/rona_pc_matrix', row.names = NULL)[,2:101]
matrix_gf_conf <- read.table('./results/offset_vs_fitness_variation/lcwp/gf_pc_matrix', row.names = NULL)[,2:101]
matrix_fitness <- read.table('./results/offset_vs_fitness_variation/lcwp/logfitness_matrix', row.names = NULL)[,2:101]


rda_conf_r2 <- c()
rona_conf_r2 <- c()
gf_conf_r2 <- c()

for (i in seq(1,30)){
  rda_conf_r2 <- c(rda_conf_r2, summary(lm(as.numeric(matrix_fitness[i,]) ~ as.numeric(matrix_rda_conf[i,])^2))$adj.r.squared)
  rona_conf_r2 <- c(rona_conf_r2, summary(lm(as.numeric(matrix_fitness[i,]) ~ as.numeric(matrix_rona_conf[i,])^2))$adj.r.squared)
  gf_conf_r2 <- c(gf_conf_r2, summary(lm(as.numeric(matrix_fitness[i,]) ~ as.numeric(matrix_gf_conf[i,])^2))$adj.r.squared)
}



result_r2_hcwp <- cbind( rda_conf_r2,  rona_conf_r2,  gf_conf_r2)
colnames(result_r2_hcwp) <- c( 'rda_pc', 'rona_pc', 'gf_pc')
write.table(result_r2_hcwp, './results/offset_vs_fitness_variation/lcwp/r2_pc_matrix')

# LCHP
matrix_rda_conf <- read.table('./results/offset_vs_fitness_variation/lchp/rda_pc_matrix', row.names = NULL)[,2:101]
matrix_rona_conf <- read.table('./results/offset_vs_fitness_variation/lchp/rona_pc_matrix', row.names = NULL)[,2:101]
matrix_gf_conf <- read.table('./results/offset_vs_fitness_variation/lchp/gf_pc_matrix', row.names = NULL)[,2:101]
matrix_fitness <- read.table('./results/offset_vs_fitness_variation/lchp/logfitness_matrix', row.names = NULL)[,2:101]


rda_conf_r2 <- c()
rona_conf_r2 <- c()
gf_conf_r2 <- c()

for (i in seq(1,30)){
  rda_conf_r2 <- c(rda_conf_r2, summary(lm(as.numeric(matrix_fitness[i,]) ~ as.numeric(matrix_rda_conf[i,])^2))$adj.r.squared)
  rona_conf_r2 <- c(rona_conf_r2, summary(lm(as.numeric(matrix_fitness[i,]) ~ as.numeric(matrix_rona_conf[i,])^2))$adj.r.squared)
  gf_conf_r2 <- c(gf_conf_r2, summary(lm(as.numeric(matrix_fitness[i,]) ~ as.numeric(matrix_gf_conf[i,])^2))$adj.r.squared)
}



result_r2_hcwp <- cbind( rda_conf_r2,  rona_conf_r2,  gf_conf_r2)
colnames(result_r2_hcwp) <- c( 'rda_pc', 'rona_pc', 'gf_pc')
write.table(result_r2_hcwp, './results/offset_vs_fitness_variation/lchp/r2_pc_matrix')