library(ggplot2)
library(qvalue)
library(LEA)
library(gradientForest)
library(dplyr)
library(RColorBrewer)

setwd('~/Documents/Rwork/diversityoffset/')
set.seed(1)

####################
# I - Collect data #
####################

# Genome data
genome_mil <- read.table("./data/data_mil/FreqAll_retainedSNP_173pop_AMMA2050_MAF10_NoNA.txt", header=TRUE, 
                         sep="\t", fileEncoding="iso-8859-1")
# Climate data
climate_mil<- read.table("./data/data_mil/ClimateMetrics_AMMApop_RealObservation.txt", 
                         header=TRUE, sep="\t", fileEncoding="iso-8859-1")

genome_align <- genome_mil[,climate_mil$Accession]
climate_only_mil <- climate_mil[,6:179]
colsd <- apply(climate_only_mil, 2, sd)
climate_only_mil <- climate_only_mil[, colsd != 0]

# We separate data from precipitation / temperature / radiation..
climate_precipitation <- climate_only_mil[,1:66]
climate_temperature <- climate_only_mil[, 67:150]
climate_radiation <- climate_only_mil[, 151:168]

colname_prec <- colnames(climate_precipitation)
colname_temp <- colnames(climate_temperature)
colname_rad <- colnames(climate_radiation)
#.. and create PC from it
pc_prec <- prcomp(climate_precipitation)
pc_temp <- prcomp(climate_temperature)
pc_rad <- prcomp(climate_radiation)
genome_align <- as.matrix(t(genome_align))

Y <- genome_align

#On doit maintenant sélectionner nos variables, on va prendre 15 axes températures
# 10 axes précipitaions et 2 axes radiation

X <- cbind(pc_prec$x[,1:10], pc_temp$x[,1:15], pc_rad$x[,1:2])



#########################
# III - Generate XPRED   #
#########################

PixelID.Sadore <- "39_15"
Pixels.value <- read.table("./data/data_mil/ClimateMetrics_pixelsINarea_RealObservation.txt",
                           header=TRUE, sep="\t", 
                           fileEncoding="iso-8859-1")

sadore.env <- Pixels.value[Pixels.value$PixelID==PixelID.Sadore,][4:177]
sadore.env <- sadore.env[, colsd != 0]
sadore.prec <- sadore.env[,1:66]
sadore.temp <- sadore.env[, 67:150]
sadore.rad <- sadore.env[, 151:168]


X.prec <- predict(pc_prec, sadore.prec)
X.temp <- predict(pc_temp, sadore.temp)
X.rad <- predict(pc_rad, sadore.rad)
X_sadore <- c(X.prec[,1:10], X.temp[,1:15], X.rad[,1:2])
X.pred <- matrix(rep(X_sadore, 173), ncol=27, byrow=T)
colnames(X.pred) <- colnames(X)

m.x <- apply(X, 2, mean)
sd.x <- apply(X, 2, sd)

X <- t(t(X) - m.x) %*% diag(1/sd.x)
X.pred <- t(t(X.pred) - m.x) %*% diag(1/sd.x)

#########################
# II - Select SNPs set #
#########################

lfmm.obj <- lfmm2(Y, X, 10, effect.sizes = T)
pv.obj <- lfmm2.test(lfmm.obj, Y, X, full=T)
qv.obj <- qvalue(pv.obj$pvalues, fdr.level=0.02)
snps.set <- which(qv.obj$significant)

#########################
# IV - Genetic gap      #
#########################

genetic_gap_scale <- genetic.gap(Y, X, pred.env=X.pred, K=10, scale=T, candidate.loci=snps.set)

barplot(genetic_gap_scale$eigenvalues, col="steelblue3", cex.axis=2)
# We multiply each vectors by its associated eigenvalues

importance <- genetic_gap_scale$vectors^2

imp_prec_1 <- sum(importance[1:10,1])
imp_temp_1 <- sum(importance[11:25,1])
imp_rad_1 <- sum(importance[26:27,1])

imp_prec_2 <- sum(importance[1:10,2])
imp_temp_2 <- sum(importance[11:25,2])
imp_rad_2 <- sum(importance[26:27,2])

imp_prec_3 <- sum(importance[1:10,3])
imp_temp_3 <- sum(importance[11:25,3])
imp_rad_3 <- sum(importance[26:27,3])

eigenvalue1 <- genetic_gap_scale$eigenvalues[1]
eigenvalue2 <- genetic_gap_scale$eigenvalues[2]
eigenvalue3 <- genetic_gap_scale$eigenvalues[3]

contrib.vec1 <- eigenvalue1 / (eigenvalue1 + eigenvalue2 + eigenvalue3)
contrib.vec2 <- eigenvalue2 / (eigenvalue1 + eigenvalue2 + eigenvalue3)
contrib.vec3 <- eigenvalue3 / (eigenvalue1 + eigenvalue2 + eigenvalue3)


imp_prec <- contrib.vec1 * imp_prec_1 + contrib.vec2 * imp_prec_2 + contrib.vec3 * imp_prec_3
imp_temp <- contrib.vec1 * imp_temp_1 + contrib.vec2 * imp_temp_2 + contrib.vec3 * imp_temp_3
imp_rad <- contrib.vec1 * imp_rad_1 + contrib.vec2 * imp_rad_2 + contrib.vec3 * imp_rad_3

df_imp <- cbind(c(imp_temp_1, imp_prec_1, imp_rad_1))
colnames(df_imp) <- c("Importance")
df_imp <- as.data.frame(df_imp)
df_imp[,"Variable"] <- c("Temperature", "Precipitation", "Radiation")

coul <- brewer.pal(3, "Set2") 
barplot(height=df_imp$Importance, names=df_imp$Variable, col=coul)



df_imp <- cbind(c(imp_temp_2, imp_prec_2, imp_rad_2))
colnames(df_imp) <- c("Importance")
df_imp <- as.data.frame(df_imp)
df_imp[,"Variable"] <- c("Temperature", "Precipitation", "Radiation")

coul <- brewer.pal(3, "Set2") 
barplot(height=df_imp$Importance, names=df_imp$Variable, col=coul)

df_imp <- cbind(c(imp_temp_3, imp_prec_3, imp_rad_3))
colnames(df_imp) <- c("Importance")
df_imp <- as.data.frame(df_imp)
df_imp[,"Variable"] <- c("Temperature", "Precipitation", "Radiation")

coul <- brewer.pal(3, "Set2") 
barplot(height=df_imp$Importance, names=df_imp$Variable, col=coul)

legend("topright", 
       legend = c("Temperature", "Precipitation", "Radiation"), 
       fill = coul)

