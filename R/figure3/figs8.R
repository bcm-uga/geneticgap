library(ggplot2)
# script to generate data for Figure S7

setwd('~/Documents/Rwork/diversityoffset/')
set.seed(1)

# HCWP
matrix_gg_hcwp <- read.table('./results/offset_vs_fitness_variation/hcwp/gg_matrix_wg', row.names = NULL)[,2:101]
matrix_rda_hcwp <- read.table('./results/offset_vs_fitness_variation/hcwp/rda_conf_matrix_wg', row.names = NULL)[,2:101]
matrix_rona_hcwp <- read.table('./results/offset_vs_fitness_variation/hcwp/rona_conf_matrix_wg', row.names = NULL)[,2:101]
matrix_gf_hcwp <- read.table('./results/offset_vs_fitness_variation/hcwp/gf_conf_matrix_wg', row.names = NULL)[,2:101]


# LCWP
matrix_gg_lcwp <- read.table('./results/offset_vs_fitness_variation/lcwp/gg_matrix_wg', row.names = NULL)[,2:101]
matrix_rda_lcwp <- read.table('./results/offset_vs_fitness_variation/lcwp/rda_conf_matrix_wg', row.names = NULL)[,2:101]
matrix_rona_lcwp <- read.table('./results/offset_vs_fitness_variation/lcwp/rona_conf_matrix_wg', row.names = NULL)[,2:101]
matrix_gf_lcwp <- read.table('./results/offset_vs_fitness_variation/lcwp/gf_conf_matrix_wg', row.names = NULL)[,2:101]


# HCHP
matrix_gg_hchp <- read.table('./results/offset_vs_fitness_variation/hchp/gg_matrix_wg', row.names = NULL)[,2:101]
matrix_rda_hchp <- read.table('./results/offset_vs_fitness_variation/hchp/rda_conf_matrix_wg', row.names = NULL)[,2:101]
matrix_rona_hchp <- read.table('./results/offset_vs_fitness_variation/hchp/rona_conf_matrix_wg', row.names = NULL)[,2:101]
matrix_gf_hchp <- read.table('./results/offset_vs_fitness_variation/hchp/gf_conf_matrix_wg', row.names = NULL)[,2:101]


# LCHP
matrix_gg_lchp <- read.table('./results/offset_vs_fitness_variation/lchp/gg_matrix_wg', row.names = NULL)[,2:101]
matrix_rda_lchp <- read.table('./results/offset_vs_fitness_variation/lchp/rda_conf_matrix_wg', row.names = NULL)[,2:101]
matrix_rona_lchp <- read.table('./results/offset_vs_fitness_variation/lchp/rona_conf_matrix_wg', row.names = NULL)[,2:101]
matrix_gf_lchp <- read.table('./results/offset_vs_fitness_variation/lchp/gf_conf_matrix_wg', row.names = NULL)[,2:101]


matrix_gg <- rbind(matrix_gg_lcwp, matrix_gg_lchp, matrix_gg_hcwp, matrix_gg_hchp)
matrix_rda <- rbind(matrix_rda_lcwp, matrix_rda_lchp, matrix_rda_hcwp, matrix_rda_hchp)
matrix_rona <- rbind(matrix_rona_lcwp, matrix_rona_lchp, matrix_rona_hcwp, matrix_rona_hchp)
matrix_gf <- rbind(matrix_gf_lcwp, matrix_gf_lchp, matrix_gf_hcwp, matrix_gf_hchp)


data <- data.frame(cbind(unlist(matrix_gg),unlist(matrix_rda)))
colnames(data) <- c("GG", "dist")

ggplot(data,aes(GG, dist)) +
  geom_point(shape=21, colour="black", fill="pink3", size=13, stroke=3) +
  geom_smooth(method='lm', se=F, col='orange1', size=6, linetype="dashed") +
  ylab("") +
  xlab("") +
  theme_bw() +
  theme(legend.position="none", 
        text = element_text(size=60),
        axis.text.x=element_blank(),
        axis.text.y=element_blank())


data <- data.frame(cbind(unlist(matrix_gg),unlist(matrix_rona)^2))
colnames(data) <- c("GG", "dist")

ggplot(data,aes(GG, dist)) +
  geom_point(shape=21, colour="black", fill="lightgoldenrod2", size=13, stroke=3) +
  geom_smooth(method='lm', se=F, col='orange1', size=6, linetype="dashed") +
  ylab("") +
  xlab("") +
  theme_bw() +
  theme(legend.position="none", 
        text = element_text(size=60),
        axis.text.x=element_blank(),
        axis.text.y=element_blank())

data <- data.frame(cbind(unlist(matrix_gg),unlist(matrix_gf)^2))
colnames(data) <- c("GG", "dist")

ggplot(data,aes(GG, dist)) +
  geom_point(shape=21, colour="black", fill="palegreen3", size=13, stroke=3) +
  geom_smooth(method='lm', se=F, col='orange1', size=6, linetype="dashed") +
  ylab("") +
  xlab("") +
  theme_bw() +
  theme(legend.position="none", 
        text = element_text(size=60),
        axis.text.x=element_blank(),
        axis.text.y=element_blank())


