library(ggplot2)
library(cowplot)

setwd('~/Documents/Rwork/diversityoffset/')

# LCWP
matrix_gg <- read.table('./results/offset_vs_fitness_variation/lcwp/gg_matrix', row.names = NULL)[,2:101]
matrix_rda <- read.table('./results/offset_vs_fitness_variation/lcwp/rda_conf_matrix', row.names = NULL)[,2:101]
matrix_rona <- read.table('./results/offset_vs_fitness_variation/lcwp/rona_conf_matrix', row.names = NULL)[,2:101]
matrix_gf <- read.table('./results/offset_vs_fitness_variation/lcwp/gf_conf_matrix', row.names = NULL)[,2:101]
matrix_dist <- read.table('./results/offset_vs_fitness_variation/lcwp/causaldist_matrix', row.names = NULL)[,2:101]
matrix_fitness <- read.table('./results/offset_vs_fitness_variation/lcwp/logfitness_matrix', row.names = NULL)[,2:101]

#SCENARIO 11

#GG
data_fitness <- data.frame(cbind(as.numeric(matrix_fitness[11,]),as.numeric(matrix_gg[11,])))
colnames(data_fitness) <- c("Fitness", "GG")

ggplot(data_fitness,aes(GG, Fitness)) +
  geom_point(shape=21, colour="black", fill="navyblue", size=13, stroke=3) +
  geom_smooth(method='lm', se=F, col='orange1', size=6, linetype="dashed") +
  ylab("") +
  xlab("") +
  theme_bw() +
  theme(legend.position="none", 
        text = element_text(size=60),
        axis.text.x=element_blank(),
        axis.text.y=element_blank())

summary(lm(as.numeric(matrix_fitness[11,]) ~ as.numeric(matrix_gg[11,])))
#RDA
data_fitness <- data.frame(cbind(as.numeric(matrix_fitness[11,]),as.numeric(matrix_rda[11,])))
colnames(data_fitness) <- c("Fitness", "GG")

ggplot(data_fitness,aes(GG, Fitness)) +
  geom_point(shape=21, colour="black", fill="pink3", size=13, stroke=3) +
  geom_smooth(method='lm', se=F, col='orange1', size=6, linetype="dashed") +
  ylab("") +
  xlab("") +
  theme_bw() +
  theme(legend.position="none", 
        text = element_text(size=60))

rda.squared <- as.numeric(matrix_rda[11,])^2
summary(lm(as.numeric(matrix_fitness[11,]) ~ rda.squared))


#GF
data_fitness <- data.frame(cbind(as.numeric(matrix_fitness[11,]),as.numeric(matrix_gf[11,])^2))
colnames(data_fitness) <- c("Fitness", "GG")

ggplot(data_fitness,aes(GG, Fitness)) +
  geom_point(shape=21, colour="black", fill="palegreen3", size=13, stroke=3) +
  geom_smooth(method='lm', se=F, col='orange1', size=6, linetype="dashed") +
  ylab("") +
  xlab("") +
  theme_bw() +
  theme(legend.position="none", 
        text = element_text(size=60))

gf.squared <- as.numeric(matrix_gf[11,])^2
summary(lm(as.numeric(matrix_fitness[11,]) ~ gf.squared))

#RONA
data_fitness <- data.frame(cbind(as.numeric(matrix_fitness[11,]),as.numeric(matrix_rona[11,])^2))
colnames(data_fitness) <- c("Fitness", "GG")

ggplot(data_fitness,aes(GG, Fitness)) +
  geom_point(shape=21, colour="black", fill="lightgoldenrod2", size=13, stroke=3) +
  geom_smooth(method='lm', se=F, col='orange1', size=6, linetype="dashed") +
  ylab("") +
  xlab("") +
  theme_bw() +
  xlim(0, 0.13) +
  theme(legend.position="none", 
        text = element_text(size=60),
        axis.text.x=element_blank(),
        axis.text.y=element_blank())

rona.squared <- as.numeric(matrix_rona[11,])^2
summary(lm(as.numeric(matrix_fitness[11,]) ~ rona.squared))