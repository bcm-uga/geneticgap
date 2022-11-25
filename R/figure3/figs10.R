library(ggplot2)
library(qvalue)
library(LEA)

setwd('~/Documents/Rwork/diversityoffset/')
set.seed(1)


# Here we need to compare genetic gap to glm
# The idea is to prove that genetic gap is a genetic distance.

data <- extract_data('./data/poly_exp/2/', MAF=0.01)

pop <- get_pop(data$pos_step1, 12, 12, 12, 12)
L <- dim(data$genome_step1)[2]

# Get snps.set

Y <- data$genome_step1
X <- cbind(data$var1_step1, data$var2_step1)
X.pred <- cbind(data$var1_pred, data$var2_pred)


haploid.obj <- haploidisation(Y, pop, X, X.pred)

Y <- haploid.obj$haploid_matrix
X <- haploid.obj$haploid_env_var
X.pred <- haploid.obj$haploid_env_var_pred
pop <- haploid.obj$haploid_population_labels
lfmm.obj <- lfmm2(Y, X, 1)
pv.obj <- lfmm2.test(lfmm.obj, Y, X, full=T)
qv.obj <- qvalue(pv.obj$pvalues, fdr.level=0.1)
snps.set <- which(qv.obj$significant)


gg.obj <- genetic.gap(Y, X, pred.env=X.pred, K=1, candidate.loci = snps.set)
gg.mean <- get_mean_by_pop(gg.obj$offset, pop)
glm.obj <- go_glm(Y, X, X.pred, pop, snps.set)
glm.obj.conf <- go_glm_conf(Y, X, X.pred, lfmm.obj@U, pop, snps.set)


# PLOT

data_glm <- data.frame(cbind(gg.mean, glm.obj.conf))
colnames(data_glm) <- c("GG", "GLM")

p2 <- ggplot(data_glm,aes(GG, GLM)) +
  geom_point(shape=21, colour="black", fill="navyblue", size=13, stroke=3) +
  geom_abline(slope=1, intercept = 0, se=F, col='indianred1', size=6, linetype="dashed")+
  ylab("") +
  xlab("") +
  theme_bw() +
  theme(legend.position="none", 
        text = element_text(size=60))

p2