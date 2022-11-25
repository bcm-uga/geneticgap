library(ggplot2)
# script to generate data for Figure S7

setwd('~/Documents/Rwork/diversityoffset/')
set.seed(1)

# HCWP
matrix_gg_hcwp <- read.table('./results/offset_vs_fitness_variation/hcwp/gg_matrix_wg', row.names = NULL)[,2:101]
matrix_dist_hcwp <- read.table('./results/offset_vs_fitness_variation/hcwp/causal_dist_matrix', row.names = NULL)[,2:101]


data <- data.frame(cbind(unlist(matrix_gg_hcwp),unlist(matrix_dist_hcwp)))
colnames(data) <- c("GG", "dist")

ggplot(data,aes(GG, dist)) +
  geom_point(shape=21, colour="black", fill="navyblue", size=13, stroke=3) +
  geom_smooth(method='lm', se=F, col='orange1', size=6, linetype="dashed") +
  ylab("") +
  xlab("") +
  theme_bw() +
  xlim(0,0.012)+
  ylim(0,0.40)+
  theme(legend.position="none", 
        text = element_text(size=60))

summary(lm(unlist(matrix_dist_hcwp) ~ unlist(matrix_gg_hcwp)))

# HCHP
matrix_gg_hchp <- read.table('./results/offset_vs_fitness_variation/hchp/gg_matrix_wg', row.names = NULL)[,2:101]
matrix_dist_hchp <- read.table('./results/offset_vs_fitness_variation/hchp/causal_dist_matrix', row.names = NULL)[,2:101]

data <- data.frame(cbind(unlist(matrix_gg_hchp),unlist(matrix_dist_hchp)))
colnames(data) <- c("GG", "dist")

ggplot(data,aes(GG, dist)) +
  geom_point(shape=21, colour="black", fill="navyblue", size=13, stroke=3) +
  geom_smooth(method='lm', se=F, col='orange1', size=6, linetype="dashed") +
  ylab("") +
  xlab("") +
  theme_bw() +
  xlim(0,0.012)+
  ylim(0,0.40)+
  theme(legend.position="none", 
        text = element_text(size=60),
        axis.text.x=element_blank(),
        axis.text.y=element_blank())

summary(lm(unlist(matrix_dist_hchp) ~ unlist(matrix_gg_hchp)))

# LCHP

matrix_gg_lchp <- read.table('./results/offset_vs_fitness_variation/lchp/gg_matrix_wg', row.names = NULL)[,2:101]
matrix_dist_lchp <- read.table('./results/offset_vs_fitness_variation/lchp/causal_dist_matrix', row.names = NULL)[,2:101]


data <- data.frame(cbind(unlist(matrix_gg_lchp),unlist(matrix_dist_lchp)))
colnames(data) <- c("GG", "dist")

ggplot(data,aes(GG, dist)) +
  geom_point(shape=21, colour="black", fill="navyblue", size=13, stroke=3) +
  geom_smooth(method='lm', se=F, col='orange1', size=6, linetype="dashed") +
  ylab("") +
  xlab("") +
  theme_bw() +
  xlim(0,0.012)+
  ylim(0,0.40)+
  theme(legend.position="none", 
        text = element_text(size=60),
        axis.text.x=element_blank(),
        axis.text.y=element_blank())

summary(lm(unlist(matrix_dist_lchp) ~ unlist(matrix_gg_lchp)))


# LCWP

matrix_gg_lcwp <- read.table('./results/offset_vs_fitness_variation/lcwp/gg_matrix_wg', row.names = NULL)[,2:101]
matrix_dist_lcwp <- read.table('./results/offset_vs_fitness_variation/lcwp/causal_dist_matrix', row.names = NULL)[,2:101]


data <- data.frame(cbind(unlist(matrix_gg_lcwp),unlist(matrix_dist_lcwp)))
colnames(data) <- c("GG", "dist")

ggplot(data,aes(GG, dist)) +
  geom_point(shape=21, colour="black", fill="navyblue", size=13, stroke=3) +
  geom_smooth(method='lm', se=F, col='orange1', size=6, linetype="dashed") +
  ylab("") +
  xlab("") +
  theme_bw() +
  xlim(0,0.012)+
  ylim(0,0.40)+
  theme(legend.position="none", 
        text = element_text(size=60),
        axis.text.x=element_blank(),
        axis.text.y=element_blank())


summary(lm(unlist(matrix_dist_lcwp) ~ unlist(matrix_gg_lcwp)))


# ALL

matrix_gg_all <- rbind(matrix_gg_hchp, matrix_gg_hcwp, matrix_gg_lchp)
matrix_dist_all <- rbind(matrix_dist_hchp, matrix_dist_hcwp, matrix_dist_lchp)


data <- data.frame(cbind(unlist(matrix_gg_all),unlist(matrix_dist_all)))
colnames(data) <- c("GG", "dist")

ggplot(data,aes(GG, dist)) +
  geom_point(shape=21, colour="black", fill="navyblue", size=13, stroke=3) +
  geom_smooth(method='lm', se=F, col='orange1', size=6, linetype="dashed") +
  ylab("") +
  xlab("") +
  theme_bw() +
  theme(legend.position="none", 
        text = element_text(size=60),
        axis.text.x=element_blank(),
        axis.text.y=element_blank())