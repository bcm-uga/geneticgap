library(ggplot2)
library(cowplot)


###################
#  no CONFOUNDING   #
###################

# HCHP
setwd('~/Documents/Rwork/diversityoffset/')
r2_matrix <- read.table('./results/offset_vs_fitness_variation/hchp/r2_matrix_wg')[,c('gg', 'rda', 'rona', 'gf')]

nb_simu <- nrow(r2_matrix)
nb_methods <- ncol(r2_matrix)

r2 <- c()
for (i in seq(1, nb_methods)){
  r2 <- c(r2, r2_matrix[,i])
}

methods <- c(rep("Genetic Gap", nb_simu),  rep("RDA", nb_simu),  rep("RONA", nb_simu),  rep("GF", nb_simu)  )

plot_df <- data.frame(cbind(methods))

colnames(plot_df) <- c("Methods")
plot_df[,"R2"] <- r2
plot_df$Methods <- factor(plot_df$Methods , levels=c("Genetic Gap", "RDA", "RONA", "GF"))

p1 <- ggplot(plot_df, aes(x=Methods, y=R2, fill=Methods)) +
  geom_boxplot() +
  xlab("") +
  ylab("") +
  scale_x_discrete("", labels=c('Genetic \n Gap',  "RDA",  "RONA", "GF")) +
  scale_fill_manual(values=c("navyblue", "pink3",  "lightgoldenrod2", "palegreen3")) +
  theme_bw() +
  ylim(0.1,0.9) +
  theme(legend.position="none", text = element_text(size=20),
        axis.text.x=element_blank(),
        axis.text.y=element_blank()) 

p1


# HCWP
setwd('~/Documents/Rwork/diversityoffset/')
r2_matrix <- read.table('./results/offset_vs_fitness_variation/hcwp/r2_matrix_wg')[,c('gg', 'rda', 'rona', 'gf')]

nb_simu <- nrow(r2_matrix)
nb_methods <- ncol(r2_matrix)

r2 <- c()
for (i in seq(1, nb_methods)){
  r2 <- c(r2, r2_matrix[,i])
}

methods <- c(rep("Genetic Gap", nb_simu),  rep("RDA", nb_simu),  rep("RONA", nb_simu),  rep("GF", nb_simu)  )

plot_df <- data.frame(cbind(methods))

colnames(plot_df) <- c("Methods")
plot_df[,"R2"] <- r2
plot_df$Methods <- factor(plot_df$Methods , levels=c("Genetic Gap", "RDA", "RONA", "GF"))

p1 <- ggplot(plot_df, aes(x=Methods, y=R2, fill=Methods)) +
  geom_boxplot() +
  xlab("") +
  ylab("") +
  scale_x_discrete("", labels=c('Genetic \n Gap',  "RDA", "RONA", "GF")) +
  scale_fill_manual(values=c("navyblue", "pink3",  "lightgoldenrod2", "palegreen3")) +
  theme_bw() +
  ylim(0.1,0.9) +
  theme(legend.position="none", text = element_text(size=20),
        axis.text.x=element_blank(),
        axis.text.y=element_blank()) 

p1



# LCHP
setwd('~/Documents/Rwork/diversityoffset/')
r2_matrix <- read.table('./results/offset_vs_fitness_variation/lchp/r2_matrix_wg')[,c('gg', 'rda', 'rona', 'gf')]

nb_simu <- nrow(r2_matrix)
nb_methods <- ncol(r2_matrix)

r2 <- c()
for (i in seq(1, nb_methods)){
  r2 <- c(r2, r2_matrix[,i])
}

methods <- c(rep("Genetic Gap", nb_simu),  rep("RDA", nb_simu),  rep("RONA", nb_simu),  rep("GF", nb_simu)  )

plot_df <- data.frame(cbind(methods))

colnames(plot_df) <- c("Methods")
plot_df[,"R2"] <- r2
plot_df$Methods <- factor(plot_df$Methods , levels=c("Genetic Gap", "RDA", "RONA", "GF" ))

p1 <- ggplot(plot_df, aes(x=Methods, y=R2, fill=Methods)) +
  geom_boxplot() +
  xlab("") +
  ylab("") +
  scale_x_discrete("", labels=c('Genetic \n Gap',  "RDA", "RONA", "GF" )) +
  scale_fill_manual(values=c("navyblue", "pink3", "lightgoldenrod2", "palegreen3")) +
  theme_bw() +
  ylim(0.1,0.9) +
  theme(legend.position="none", text = element_text(size=20),
        axis.text.x=element_blank(),
        axis.text.y=element_blank()) 

p1


# LCWP
setwd('~/Documents/Rwork/diversityoffset/')
r2_matrix <- read.table('./results/offset_vs_fitness_variation/lcwp/r2_matrix_wg')[,c('gg', 'rda', 'rona', 'gf')]

nb_simu <- nrow(r2_matrix)
nb_methods <- ncol(r2_matrix)

r2 <- c()
for (i in seq(1, nb_methods)){
  r2 <- c(r2, r2_matrix[,i])
}

methods <- c(rep("Genetic Gap", nb_simu),  rep("RDA", nb_simu),  rep("RONA", nb_simu),  rep("GF", nb_simu)  )

plot_df <- data.frame(cbind(methods))

colnames(plot_df) <- c("Methods")
plot_df[,"R2"] <- r2
plot_df$Methods <- factor(plot_df$Methods , levels=c("Genetic Gap", "RDA", "RONA", "GF" ))

p1 <- ggplot(plot_df, aes(x=Methods, y=R2, fill=Methods)) +
  geom_boxplot() +
  xlab("") +
  ylab("") +
  scale_x_discrete("", labels=c('Genetic \n Gap',  "RDA", "RONA", "GF" )) +
  scale_fill_manual(values=c("navyblue", "pink3", "lightgoldenrod2", "palegreen3")) +
  theme_bw() +
  ylim(0.1,0.9) +
  theme(legend.position="none", text = element_text(size=20),
        axis.text.x=element_blank(),
        axis.text.y=element_blank()) 

p1


###################
#   CONFOUNDING   #
###################

# HCHP
setwd('~/Documents/Rwork/diversityoffset/')
r2_matrix <- read.table('./results/offset_vs_fitness_variation/hchp/r2_matrix_wg')[,c('gg', 'rda_conf', 'rona_conf', 'gf_conf')]

nb_simu <- nrow(r2_matrix)
nb_methods <- ncol(r2_matrix)

r2 <- c()
for (i in seq(1, nb_methods)){
  r2 <- c(r2, r2_matrix[,i])
}

methods <- c(rep("Genetic Gap", nb_simu),  rep("RDA", nb_simu),  rep("RONA", nb_simu),  rep("GF", nb_simu)  )

plot_df <- data.frame(cbind(methods))

colnames(plot_df) <- c("Methods")
plot_df[,"R2"] <- r2
plot_df$Methods <- factor(plot_df$Methods , levels=c("Genetic Gap", "RDA", "RONA", "GF"))

p1 <- ggplot(plot_df, aes(x=Methods, y=R2, fill=Methods)) +
  geom_boxplot() +
  xlab("") +
  ylab("") +
  scale_x_discrete("", labels=c('Genetic \n Gap',  "RDA",  "RONA", "GF")) +
  scale_fill_manual(values=c("navyblue", "pink3",  "lightgoldenrod2", "palegreen3")) +
  theme_bw() +
  ylim(0.1,0.9) +
  theme(legend.position="none", text = element_text(size=20),
        axis.text.x=element_blank(),
        axis.text.y=element_blank()) 

p1


# HCWP
setwd('~/Documents/Rwork/diversityoffset/')
r2_matrix <- read.table('./results/offset_vs_fitness_variation/hcwp/r2_matrix_wg')[,c('gg', 'rda_conf', 'rona_conf', 'gf_conf')]

nb_simu <- nrow(r2_matrix)
nb_methods <- ncol(r2_matrix)

r2 <- c()
for (i in seq(1, nb_methods)){
  r2 <- c(r2, r2_matrix[,i])
}

methods <- c(rep("Genetic Gap", nb_simu),  rep("RDA", nb_simu),  rep("RONA", nb_simu),  rep("GF", nb_simu)  )

plot_df <- data.frame(cbind(methods))

colnames(plot_df) <- c("Methods")
plot_df[,"R2"] <- r2
plot_df$Methods <- factor(plot_df$Methods , levels=c("Genetic Gap", "RDA", "RONA", "GF"))

p1 <- ggplot(plot_df, aes(x=Methods, y=R2, fill=Methods)) +
  geom_boxplot() +
  xlab("") +
  ylab("") +
  scale_x_discrete("", labels=c('Genetic \n Gap',  "RDA", "RONA", "GF")) +
  scale_fill_manual(values=c("navyblue", "pink3",  "lightgoldenrod2", "palegreen3")) +
  theme_bw() +
  ylim(0.1,0.9) +
  theme(legend.position="none", text = element_text(size=20),
        axis.text.x=element_blank(),
        axis.text.y=element_blank()) 

p1



# LCHP
setwd('~/Documents/Rwork/diversityoffset/')
r2_matrix <- read.table('./results/offset_vs_fitness_variation/lchp/r2_matrix_wg')[,c('gg', 'rda_conf', 'rona_conf', 'gf_conf')]

nb_simu <- nrow(r2_matrix)
nb_methods <- ncol(r2_matrix)

r2 <- c()
for (i in seq(1, nb_methods)){
  r2 <- c(r2, r2_matrix[,i])
}

methods <- c(rep("Genetic Gap", nb_simu),  rep("RDA", nb_simu),  rep("RONA", nb_simu),  rep("GF", nb_simu)  )

plot_df <- data.frame(cbind(methods))

colnames(plot_df) <- c("Methods")
plot_df[,"R2"] <- r2
plot_df$Methods <- factor(plot_df$Methods , levels=c("Genetic Gap", "RDA", "RONA", "GF" ))

p1 <- ggplot(plot_df, aes(x=Methods, y=R2, fill=Methods)) +
  geom_boxplot() +
  xlab("") +
  ylab("") +
  scale_x_discrete("", labels=c('Genetic \n Gap',  "RDA", "RONA", "GF" )) +
  scale_fill_manual(values=c("navyblue", "pink3", "lightgoldenrod2", "palegreen3")) +
  theme_bw() +
  ylim(0.1,0.9) +
  theme(legend.position="none", text = element_text(size=20),
        axis.text.x=element_blank(),
        axis.text.y=element_blank()) 

p1


# LCWP
setwd('~/Documents/Rwork/diversityoffset/')
r2_matrix <- read.table('./results/offset_vs_fitness_variation/lcwp/r2_matrix_wg')[,c('gg', 'rda_conf', 'rona_conf', 'gf_conf')]

nb_simu <- nrow(r2_matrix)
nb_methods <- ncol(r2_matrix)

r2 <- c()
for (i in seq(1, nb_methods)){
  r2 <- c(r2, r2_matrix[,i])
}

methods <- c(rep("Genetic Gap", nb_simu),  rep("RDA", nb_simu),  rep("RONA", nb_simu),  rep("GF", nb_simu)  )

plot_df <- data.frame(cbind(methods))

colnames(plot_df) <- c("Methods")
plot_df[,"R2"] <- r2
plot_df$Methods <- factor(plot_df$Methods , levels=c("Genetic Gap", "RDA", "RONA", "GF" ))

p1 <- ggplot(plot_df, aes(x=Methods, y=R2, fill=Methods)) +
  geom_boxplot() +
  xlab("") +
  ylab("") +
  scale_x_discrete("", labels=c('Genetic \n Gap',  "RDA", "RONA", "GF" )) +
  scale_fill_manual(values=c("navyblue", "pink3", "lightgoldenrod2", "palegreen3")) +
  theme_bw() +
  ylim(0.1,0.9) +
  theme(legend.position="none", text = element_text(size=20),
        axis.text.x=element_blank(),
        axis.text.y=element_blank()) 

p1


