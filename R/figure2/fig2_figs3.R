# We will generate fig 2 and figS2
# In this example we will use scenario 10 from 'north migration' scenario


library(ggplot2)
library(fields)
library(qvalue)
library(lmtest)
library(LEA)

set.seed(5)

setwd("~/Documents/Rwork/diversityoffset/R/figure2/")
data <- extract_data('../../data/poly_small/1/', MAF=0.01)
pop <- get_pop(data$pos_step1, 12,12,12,12)
n <- nrow(data$genome_step1)


################
#     Fig2C    #
################

pc.obj <- prcomp(data$genome_step1)
df.plot <- as.data.frame(cbind(pc.obj$x[,1], pc.obj$x[,2],data$var1_step1, data$var2_step1))
colnames(df.plot) <- c("PC1", "PC2", "VAR1","VAR2")


png("fig2C_var2.png")
ggplot(df.plot, aes(x = PC1, y = PC2, colour = VAR2)) +
  geom_point() +
  scale_colour_gradientn(colours = terrain.colors(10)) +
  theme_bw() +
  xlab("") +
  ylab("") +
  theme( text = element_text(size=25))
dev.off()

png("fig2C_var1.png")
ggplot(df.plot, aes(x = PC1, y = PC2, colour = VAR1)) +
  geom_point() +
  scale_colour_gradientn(colours = terrain.colors(10)) +
  theme_bw() +
  xlab("") +
  ylab("") +
  theme( text = element_text(size=25)) 
dev.off()

# Correlation between PC1 and environmental variables

print("Correlation between PC1 and VAR1 :")
print(cor(pc.obj$x[,1], data$var1_step1))

print("Correlation between PC2 and VAR1 :")
print(cor(pc.obj$x[,2], data$var1_step1))

print("Correlation between PC1 and VAR2 :")
print(cor(pc.obj$x[,1], data$var2_step1))

print("Correlation between PC2 and VAR2 :")
print(cor(pc.obj$x[,2], data$var2_step1))


################
#     Fig2B    #
################

# First we will generate 2 new variables that are linear combinations of the first two predictors
# We will also create their future state (pred)
data$var3_step1 <- 0.5 * (data$var1_step1 + data$var2_step1) + rnorm(n, 0, 0.05)
data$var4_step1 <- 0.8 * data$var1_step1  - 0.3 * data$var2_step1 + 0.3 + rnorm(n, 0, 0.05)

# We will generate vector of noise for future variables

noise_var3 <- rep(1, n)
noise_var4 <- rep(1, n)
for (i in seq(1,100)){
  nb_pop <- sum(pop==i)
  noise_var3[pop==i] = rep(runif(1,-0.3, 0.3), nb_pop)
  noise_var4[pop==i] = rep(runif(1,-0.3, 0.3), nb_pop)
}
data$var3_pred <- data$var3_step1 + noise_var3
data$var4_pred <- data$var4_step1 + noise_var4


# Compute correlation between environmental variables

print(" Correlation between VAR3 and VAR1 :")
print(cor(data$var3_step1, data$var1_step1))

print(" Correlation between VAR3 and VAR2 :")
print(cor(data$var3_step1, data$var2_step1))

print(" Correlation between VAR4 and VAR1 :")
print(cor(data$var4_step1, data$var1_step1))

print(" Correlation between VAR4 and VAR2 :")
print(cor(data$var4_step1, data$var2_step1))


# We can now generate all the maps
coord <- cbind(data$pos_step1$x, data$pos_step1$y)

fit = Krig(coord, data$var1_step1, m = 2, theta = 10)
surface(fit, extrap = TRUE, main = "", cex.axis=2, levels = c(3), col=turbo(256), xlab="", ylab="", zlim=c(-0.2,1.2), xaxt="n", yaxt="n")

fit = Krig(coord, data$var2_step1, m = 2, theta = 10)
surface(fit, extrap = TRUE, main = "", cex=2, levels = c(3), col=turbo(256), xlab="", ylab="", zlim=c(-0.2,1.2), xaxt="n", yaxt="n")

fit = Krig(coord, data$var3_step1, m = 2, theta = 10)
surface(fit, extrap = TRUE, main = "", cex=2, levels = c(3), col=turbo(256), xlab="", ylab="", zlim=c(-0.2,1.2) , xaxt="n", yaxt="n")

fit = Krig(coord, data$var4_step1, m = 2, theta = 10)
surface(fit, extrap = TRUE, main = "", cex=2, levels = c(3), col=turbo(256), xlab="", ylab="", zlim=c(-0.2,1.2), xaxt="n", yaxt="n")


fit = Krig(coord, data$var1_pred, m = 2, theta = 10)
surface(fit, extrap = TRUE, main = "", cex=2, levels = c(3), col=turbo(256), xlab="", ylab="", zlim=c(-0.2,1.2), xaxt="n", yaxt="n")

fit = Krig(coord, data$var2_pred, m = 2, theta = 10)
surface(fit, extrap = TRUE, main = "", cex=2, levels = c(3), col=turbo(256), xlab="", ylab="", zlim=c(-0.2,1.2), xaxt="n", yaxt="n")

fit = Krig(coord, data$var3_pred, m = 2, theta = 10)
surface(fit, extrap = TRUE, main = "", cex=2, levels = c(3), col=turbo(256), xlab="", ylab="", zlim=c(-0.2,1.2), xaxt="n", yaxt="n")

fit = Krig(coord, data$var4_pred, m = 2, theta = 10)
surface(fit, extrap = TRUE, main = "", cex=2, levels = c(3), col=turbo(256), xlab="", ylab="", zlim=c(-0.2,1.2), xaxt="n", yaxt="n")



################
#     Fig2E    #
################

# Focus on genetic gap on the whole genome

# Select SNPs on which to evaluate the method
Y <- data$genome_step1
X <- cbind(data$var1_step1, data$var2_step1, data$var3_step1, data$var4_step1)
X.pred <- cbind(data$var1_pred, data$var2_pred, data$var3_pred, data$var4_pred)
lfmm.obj <- lfmm2(Y, X, 10)
pv.obj <- lfmm2.test(lfmm.obj, Y, X, T)
qv.obj <- qvalue(pv.obj$pvalues, fdr.level=0.05)
snps.set <- which(qv.obj$significant)

# Plot fitness_var vs genetic.gap

gg.obj <- genetic.gap(Y, X, pred.env = X.pred, K=10, scale=T, candidate.loci = snps.set)

gg <- get_mean_by_pop(gg.obj$offset, pop)
fitness_var <- get_fitness_variation(data, pop)
fitness_pred <- get_mean_by_pop(data$fitness_pred, pop)

data_fitness <- data.frame(cbind(-log(fitness_pred), gg))
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

lm_gg <- lm(-log(fitness_pred) ~ gg)
summary(lm(-log(fitness_pred) ~ gg))



causal_dist <- causal_distance(X[,1:2], X.pred[,1:2])
mean_causal_dist <- get_mean_by_pop(causal_dist, pop)

euclid_dist <- euclidean_distance(X, X.pred)
mean_euclid_dist <- get_mean_by_pop(euclid_dist, pop)


data_dist <- data.frame(cbind(mean_causal_dist, gg))
colnames(data_dist) <- c("Distance", "GG")

ggplot(data_dist,aes(GG, Distance)) +
  geom_point(shape=21, colour="black", fill="navyblue", size=13, stroke=3) +
  geom_smooth(method='lm', se=F, col='orange1', size=6, linetype="dashed") +
  ylab("") +
  xlab("") +
  theme_bw() +
  theme(legend.position="none", 
        text = element_text(size=60))

summary(lm(mean_causal_dist ~ gg))

lm_eucli <- lm(-log(fitness_pred) ~ mean_euclid_dist)
summary(lm(-log(fitness_pred) ~ mean_euclid_dist))

jtest(lm_gg, lm_eucli)
# Plot screeplot

screeplot_values <- gg.obj$eigenvalues
vector_values <- gg.obj$vectors

barplot(screeplot_values, col="steelblue3", cex.axis=2)

plot(screeplot_values, type='h', col="blue", axes=F, ylab="Eigenvalues")
axis(1, at = seq(1,4))
axis(2, at = seq(0, 0.4, by =0.05))
points(seq(1,4),screeplot_values, pch=19, col="red")


## FIG S2

## VS Distance euclidienne au carré

euclidean_distance <- function(var_cur, var_pred){
  
  nb_var <- nrow(var_cur)
  
  climate_distance <- c()
  for (i in seq(1,nb_var)){
    climate_distance <- c(climate_distance, dist(rbind(var_cur[i,], var_pred[i,])))
  }
  
  return(climate_distance)
}

X.cur <- cbind(get_mean_by_pop(data$var1_step1,pop), get_mean_by_pop(data$var2_step1, pop), get_mean_by_pop(data$var3_step1, pop), get_mean_by_pop(data$var4_step1,pop))
X.pred <- cbind(get_mean_by_pop(data$var1_pred,pop), get_mean_by_pop(data$var2_pred, pop), get_mean_by_pop(data$var3_pred,pop), get_mean_by_pop(data$var4_pred,pop))

dist_eucl <- euclidean_distance(X.cur, X.pred)**2

summary(lm(-log(fitness_pred) ~ dist_eucl))

data_eucl <- data.frame(cbind(-log(fitness_pred), dist_eucl))
colnames(data_eucl) <- c("Fitness", "Dist")

ggplot(data_eucl,aes(Dist, Fitness)) +
  geom_point(shape=21, colour="black", fill="navyblue", size=15, stroke=3) +
  geom_smooth(method='lm', se=F, col='yellow1', size=3) +
  ylab("") +
  xlab("") +
  theme_bw() +
  theme(legend.position="none", 
        text = element_text(size=60),
        axis.text.x=element_blank(),
        axis.text.y=element_blank())




# Distance causale
# (1/2)*((((optimum_var1 - phenotype_var1)/sigma_K)^2) + ((optimum_var2 - phenotype_var2)/sigma_K)^2)
# ou optimum_var1 = X.causal.pred et phenotype_var1 = X.causal.cur

causal_distance <- function(X.causal.cur, X.causal.pred){
  nb_var <- nrow(X.causal.cur)
  
  causal_distance <- c()
  for (i in seq(1,nb_var)){
    optimum_var1 <- X.causal.pred[i,1]
    phenotype_var1 <- X.causal.cur[i,1]
    optimum_var2 <- X.causal.pred[i,2]
    phenotype_var2 <- X.causal.cur[i,2]
    sigma_K <- 0.5
    causal_distance <- c(causal_distance, (1/2)*((((optimum_var1 - phenotype_var1)/sigma_K)^2) + ((optimum_var2 - phenotype_var2)/sigma_K)^2))
  }
  
  return(causal_distance)
}

X.causal.cur <- cbind(get_mean_by_pop(data$var1_step1,pop), get_mean_by_pop(data$var2_step1, pop))
X.causal.pred <- cbind(get_mean_by_pop(data$var1_pred,pop), get_mean_by_pop(data$var2_pred, pop))

causal_dist <- causal_distance(X.causal.cur, X.causal.pred)

data_caus <- data.frame(cbind(-log(fitness_pred), causal_dist))
colnames(data_caus) <- c("Fitness", "Dist")

ggplot(data_caus,aes(Dist, Fitness)) +
  geom_point(shape=21, colour="black", fill="navyblue", size=15, stroke=3) +
  geom_smooth(method='lm', se=F, col='yellow1', size=3) +
  ylab("") +
  xlab("") +
  theme_bw() +
  theme(legend.position="none", 
        text = element_text(size=60),
        axis.text.x=element_blank(),
        axis.text.y=element_blank())


################
#     Fig2D    #
################

fitness_var_map <- -log(data$fitness_pred)
fit = Krig(coord, fitness_var_map, m = 2, theta = 10)
surface(fit, extrap = T, main = "", cex.axis = 2, levels = c(3), col=turbo(256), xlab="", ylab="", type="C", xaxt="n", yaxt="n")

# We need one offset value per individual 



fit = Krig(coord, gg.obj$offset, m = 2, theta = 10)
surface(fit, extrap = T, main = "", cex.main = 2, levels = c(3), col=turbo(256), xlab="", ylab="", type="C",  xaxt="n", yaxt="n")
