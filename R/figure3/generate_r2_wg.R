# The aim here is to generate R2 for the different scenario


setwd('~/Documents/Rwork/diversityoffset/')
set.seed(1)


# HCWP
matrix_gg <- read.table('./results/offset_vs_fitness_variation/hcwp/gg_matrix_wg', row.names = NULL)[,2:101]
matrix_rda <- read.table('./results/offset_vs_fitness_variation/hcwp/rda_matrix_wg', row.names = NULL)[,2:101]
matrix_rda_conf <- read.table('./results/offset_vs_fitness_variation/hcwp/rda_conf_matrix_wg', row.names = NULL)[,2:101]
matrix_rona <- read.table('./results/offset_vs_fitness_variation/hcwp/rona_matrix_wg', row.names = NULL)[,2:101]
matrix_rona_conf <- read.table('./results/offset_vs_fitness_variation/hcwp/rona_conf_matrix_wg', row.names = NULL)[,2:101]
matrix_gf <- read.table('./results/offset_vs_fitness_variation/hcwp/gf_matrix_wg', row.names = NULL)[,2:101]
matrix_gf_conf <- read.table('./results/offset_vs_fitness_variation/hcwp/gf_conf_matrix_wg', row.names = NULL)[,2:101]
matrix_dist <- read.table('./results/offset_vs_fitness_variation/hcwp/causaldist_matrix_wg', row.names = NULL)[,2:101]
matrix_fitness <- read.table('./results/offset_vs_fitness_variation/hcwp/logfitness_matrix_wg', row.names = NULL)[,2:101]

gg_r2 <- c()
rda_r2 <- c()
rda_conf_r2 <- c()
rona_r2 <- c()
rona_conf_r2 <- c()
gf_r2 <- c()
gf_conf_r2 <- c()
dist_r2 <- c()

for (i in seq(1,30)){
  gg_r2 <- c(gg_r2, summary(lm(as.numeric(matrix_fitness[i,]) ~ as.numeric(matrix_gg[i,])))$adj.r.squared)
  rda_r2 <- c(rda_r2, summary(lm(as.numeric(matrix_fitness[i,]) ~ as.numeric(matrix_rda[i,])^2))$adj.r.squared)
  rda_conf_r2 <- c(rda_conf_r2, summary(lm(as.numeric(matrix_fitness[i,]) ~ as.numeric(matrix_rda_conf[i,])^2))$adj.r.squared)
  rona_r2 <- c(rona_r2, summary(lm(as.numeric(matrix_fitness[i,]) ~ as.numeric(matrix_rona[i,])^2))$adj.r.squared)
  rona_conf_r2 <- c(rona_conf_r2, summary(lm(as.numeric(matrix_fitness[i,]) ~ as.numeric(matrix_rona_conf[i,])^2))$adj.r.squared)
  gf_r2 <- c(gf_r2, summary(lm(as.numeric(matrix_fitness[i,]) ~ as.numeric(matrix_gf[i,])^2))$adj.r.squared)
  gf_conf_r2 <- c(gf_conf_r2, summary(lm(as.numeric(matrix_fitness[i,]) ~ as.numeric(matrix_gf_conf[i,])^2))$adj.r.squared)
  dist_r2 <- c(dist_r2, summary(lm(as.numeric(matrix_fitness[i,]) ~ as.numeric(matrix_dist[i,])))$adj.r.squared)
}



result_r2_hcwp <- cbind(gg_r2,  rda_r2, rda_conf_r2, rona_r2, rona_conf_r2,  gf_r2, gf_conf_r2, dist_r2)
colnames(result_r2_hcwp) <- c('gg','rda', 'rda_conf', 'rona', 'rona_conf', 'gf', 'gf_conf', 'distance')
write.table(result_r2_hcwp, './results/offset_vs_fitness_variation/hcwp/r2_matrix_wg')



# HCHP
matrix_gg <- read.table('./results/offset_vs_fitness_variation/hchp/gg_matrix_wg', row.names = NULL)[,2:101]
matrix_rda <- read.table('./results/offset_vs_fitness_variation/hchp/rda_matrix_wg', row.names = NULL)[,2:101]
matrix_rda_conf <- read.table('./results/offset_vs_fitness_variation/hchp/rda_conf_matrix_wg', row.names = NULL)[,2:101]
matrix_rona <- read.table('./results/offset_vs_fitness_variation/hchp/rona_matrix_wg', row.names = NULL)[,2:101]
matrix_rona_conf <- read.table('./results/offset_vs_fitness_variation/hchp/rona_conf_matrix_wg', row.names = NULL)[,2:101]
matrix_gf <- read.table('./results/offset_vs_fitness_variation/hchp/gf_matrix_wg', row.names = NULL)[,2:101]
matrix_gf_conf <- read.table('./results/offset_vs_fitness_variation/hchp/gf_conf_matrix_wg', row.names = NULL)[,2:101]
matrix_dist <- read.table('./results/offset_vs_fitness_variation/hchp/causaldist_matrix_wg', row.names = NULL)[,2:101]
matrix_fitness <- read.table('./results/offset_vs_fitness_variation/hchp/logfitness_matrix_wg', row.names = NULL)[,2:101]

gg_r2 <- c()
rda_r2 <- c()
rda_conf_r2 <- c()
rona_r2 <- c()
rona_conf_r2 <- c()
gf_r2 <- c()
gf_conf_r2 <- c()
dist_r2 <- c()

for (i in seq(1,30)){
  gg_r2 <- c(gg_r2, summary(lm(as.numeric(matrix_fitness[i,]) ~ as.numeric(matrix_gg[i,])))$adj.r.squared)
  rda_r2 <- c(rda_r2, summary(lm(as.numeric(matrix_fitness[i,]) ~ as.numeric(matrix_rda[i,])^2))$adj.r.squared)
  rda_conf_r2 <- c(rda_conf_r2, summary(lm(as.numeric(matrix_fitness[i,]) ~ as.numeric(matrix_rda_conf[i,])^2))$adj.r.squared)
  rona_r2 <- c(rona_r2, summary(lm(as.numeric(matrix_fitness[i,]) ~ as.numeric(matrix_rona[i,])^2))$adj.r.squared)
  rona_conf_r2 <- c(rona_conf_r2, summary(lm(as.numeric(matrix_fitness[i,]) ~ as.numeric(matrix_rona_conf[i,])^2))$adj.r.squared)
  gf_r2 <- c(gf_r2, summary(lm(as.numeric(matrix_fitness[i,]) ~ as.numeric(matrix_gf[i,])^2))$adj.r.squared)
  gf_conf_r2 <- c(gf_conf_r2, summary(lm(as.numeric(matrix_fitness[i,]) ~ as.numeric(matrix_gf_conf[i,])^2))$adj.r.squared)
  dist_r2 <- c(dist_r2, summary(lm(as.numeric(matrix_fitness[i,]) ~ as.numeric(matrix_dist[i,])))$adj.r.squared)
}



result_r2_hcwp <- cbind(gg_r2,  rda_r2, rda_conf_r2, rona_r2, rona_conf_r2,  gf_r2, gf_conf_r2, dist_r2)
colnames(result_r2_hcwp) <- c('gg','rda', 'rda_conf', 'rona', 'rona_conf', 'gf', 'gf_conf', 'distance')
write.table(result_r2_hcwp, './results/offset_vs_fitness_variation/hchp/r2_matrix_wg')



# LCWP
matrix_gg <- read.table('./results/offset_vs_fitness_variation/lcwp/gg_matrix_wg', row.names = NULL)[,2:101]
matrix_rda <- read.table('./results/offset_vs_fitness_variation/lcwp/rda_matrix_wg', row.names = NULL)[,2:101]
matrix_rda_conf <- read.table('./results/offset_vs_fitness_variation/lcwp/rda_conf_matrix_wg', row.names = NULL)[,2:101]
matrix_rona <- read.table('./results/offset_vs_fitness_variation/lcwp/rona_matrix_wg', row.names = NULL)[,2:101]
matrix_rona_conf <- read.table('./results/offset_vs_fitness_variation/lcwp/rona_conf_matrix_wg', row.names = NULL)[,2:101]
matrix_gf <- read.table('./results/offset_vs_fitness_variation/lcwp/gf_matrix_wg', row.names = NULL)[,2:101]
matrix_gf_conf <- read.table('./results/offset_vs_fitness_variation/lcwp/gf_conf_matrix_wg', row.names = NULL)[,2:101]
matrix_dist <- read.table('./results/offset_vs_fitness_variation/lcwp/causaldist_matrix_wg', row.names = NULL)[,2:101]
matrix_fitness <- read.table('./results/offset_vs_fitness_variation/lcwp/logfitness_matrix_wg', row.names = NULL)[,2:101]

gg_r2 <- c()
rda_r2 <- c()
rda_conf_r2 <- c()
rona_r2 <- c()
rona_conf_r2 <- c()
gf_r2 <- c()
gf_conf_r2 <- c()
dist_r2 <- c()

for (i in seq(1,30)){
  gg_r2 <- c(gg_r2, summary(lm(as.numeric(matrix_fitness[i,]) ~ as.numeric(matrix_gg[i,])))$adj.r.squared)
  rda_r2 <- c(rda_r2, summary(lm(as.numeric(matrix_fitness[i,]) ~ as.numeric(matrix_rda[i,])^2))$adj.r.squared)
  rda_conf_r2 <- c(rda_conf_r2, summary(lm(as.numeric(matrix_fitness[i,]) ~ as.numeric(matrix_rda_conf[i,])^2))$adj.r.squared)
  rona_r2 <- c(rona_r2, summary(lm(as.numeric(matrix_fitness[i,]) ~ as.numeric(matrix_rona[i,])^2))$adj.r.squared)
  rona_conf_r2 <- c(rona_conf_r2, summary(lm(as.numeric(matrix_fitness[i,]) ~ as.numeric(matrix_rona_conf[i,])^2))$adj.r.squared)
  gf_r2 <- c(gf_r2, summary(lm(as.numeric(matrix_fitness[i,]) ~ as.numeric(matrix_gf[i,])^2))$adj.r.squared)
  gf_conf_r2 <- c(gf_conf_r2, summary(lm(as.numeric(matrix_fitness[i,]) ~ as.numeric(matrix_gf_conf[i,])^2))$adj.r.squared)
  dist_r2 <- c(dist_r2, summary(lm(as.numeric(matrix_fitness[i,]) ~ as.numeric(matrix_dist[i,])))$adj.r.squared)
}



result_r2_hcwp <- cbind(gg_r2,  rda_r2, rda_conf_r2, rona_r2, rona_conf_r2,  gf_r2, gf_conf_r2, dist_r2)
colnames(result_r2_hcwp) <- c('gg','rda', 'rda_conf', 'rona', 'rona_conf', 'gf', 'gf_conf', 'distance')
write.table(result_r2_hcwp, './results/offset_vs_fitness_variation/lcwp/r2_matrix_wg')




# LCHP
matrix_gg <- read.table('./results/offset_vs_fitness_variation/lchp/gg_matrix_wg', row.names = NULL)[,2:101]
matrix_rda <- read.table('./results/offset_vs_fitness_variation/lchp/rda_matrix_wg', row.names = NULL)[,2:101]
matrix_rda_conf <- read.table('./results/offset_vs_fitness_variation/lchp/rda_conf_matrix_wg', row.names = NULL)[,2:101]
matrix_rona <- read.table('./results/offset_vs_fitness_variation/lchp/rona_matrix_wg', row.names = NULL)[,2:101]
matrix_rona_conf <- read.table('./results/offset_vs_fitness_variation/lchp/rona_conf_matrix_wg', row.names = NULL)[,2:101]
matrix_gf <- read.table('./results/offset_vs_fitness_variation/lchp/gf_matrix_wg', row.names = NULL)[,2:101]
matrix_gf_conf <- read.table('./results/offset_vs_fitness_variation/lchp/gf_conf_matrix_wg', row.names = NULL)[,2:101]
matrix_dist <- read.table('./results/offset_vs_fitness_variation/lchp/causaldist_matrix_wg', row.names = NULL)[,2:101]
matrix_fitness <- read.table('./results/offset_vs_fitness_variation/lchp/logfitness_matrix_wg', row.names = NULL)[,2:101]

gg_r2 <- c()
rda_r2 <- c()
rda_conf_r2 <- c()
rona_r2 <- c()
rona_conf_r2 <- c()
gf_r2 <- c()
gf_conf_r2 <- c()
dist_r2 <- c()

for (i in seq(1,30)){
  gg_r2 <- c(gg_r2, summary(lm(as.numeric(matrix_fitness[i,]) ~ as.numeric(matrix_gg[i,])))$adj.r.squared)
  rda_r2 <- c(rda_r2, summary(lm(as.numeric(matrix_fitness[i,]) ~ as.numeric(matrix_rda[i,])^2))$adj.r.squared)
  rda_conf_r2 <- c(rda_conf_r2, summary(lm(as.numeric(matrix_fitness[i,]) ~ as.numeric(matrix_rda_conf[i,])^2))$adj.r.squared)
  rona_r2 <- c(rona_r2, summary(lm(as.numeric(matrix_fitness[i,]) ~ as.numeric(matrix_rona[i,])^2))$adj.r.squared)
  rona_conf_r2 <- c(rona_conf_r2, summary(lm(as.numeric(matrix_fitness[i,]) ~ as.numeric(matrix_rona_conf[i,])^2))$adj.r.squared)
  gf_r2 <- c(gf_r2, summary(lm(as.numeric(matrix_fitness[i,]) ~ as.numeric(matrix_gf[i,])^2))$adj.r.squared)
  gf_conf_r2 <- c(gf_conf_r2, summary(lm(as.numeric(matrix_fitness[i,]) ~ as.numeric(matrix_gf_conf[i,])^2))$adj.r.squared)
  dist_r2 <- c(dist_r2, summary(lm(as.numeric(matrix_fitness[i,]) ~ as.numeric(matrix_dist[i,])))$adj.r.squared)
}



result_r2_hcwp <- cbind(gg_r2,  rda_r2, rda_conf_r2, rona_r2, rona_conf_r2,  gf_r2, gf_conf_r2, dist_r2)
colnames(result_r2_hcwp) <- c('gg','rda', 'rda_conf', 'rona', 'rona_conf', 'gf', 'gf_conf', 'distance')
write.table(result_r2_hcwp, './results/offset_vs_fitness_variation/lchp/r2_matrix_wg')




