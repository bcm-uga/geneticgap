library(ggplot2)
library(qvalue)
library(LEA)
library(dplyr)

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


# relate genetic gap to log(fitness_pred)

pheno <- read.csv('./data/data_mil/RAWDATA_Pheno_ProxyYield_2016-2017_6Trials_Sadore.csv', header=T, sep=";", na.string="NA")
PhenoMoyAccess <- pheno %>% 
  # Compute mean phenotypic measures per accession
  group_by(Accession) %>% 
  summarise(meanNCS=mean(NoCollSpike, na.rm=TRUE), sdNCS=sd(NoCollSpike, na.rm=TRUE),
            meanTWS=mean(TWeightSeed, na.rm=TRUE), sdTWS=sd(TWeightSeed, na.rm=TRUE),
            meanW100=mean(Weight100, na.rm=TRUE), sdW100=sd(Weight100, na.rm=TRUE),
            meanWSS=mean(TWeightSpike, na.rm=TRUE), sdWSS=sd(TWeightSpike, na.rm=TRUE)
  )
PhenoMoyAccess <- as.data.frame(PhenoMoyAccess)
df.offset <- data.frame(rownames(genome_align))
colnames(df.offset) <- c("Accession")
df.offset["offsetcausal"] <- genetic_gap_scale$offset
PhenoMoyAccess <- merge(PhenoMoyAccess, df.offset)
remove_outlier <- log(PhenoMoyAccess$meanTWS) > 3

# And plot

data_lfmm <- data.frame(cbind(PhenoMoyAccess$offsetcausal[remove_outlier], log(PhenoMoyAccess$meanTWS[remove_outlier])))
colnames(data_lfmm) <- c("LFMM", "Fitness")


p3 <- ggplot(data_lfmm,aes(LFMM, Fitness)) +
  geom_point(shape=21, colour="black", fill="navyblue", size=13, stroke=3) +
  geom_smooth(method='lm', se=F, col='orange1', size=6, linetype="dashed") +
  ylab("") +
  xlab("") +
  theme_bw() +
  theme(legend.position="none", 
        text = element_text(size=60),
        axis.text.x=element_blank(),
        axis.text.y=element_blank())


p3

lm_gg <- lm(log(PhenoMoyAccess$meanTWS[remove_outlier]) ~ PhenoMoyAccess$offsetcausal[remove_outlier])
summary(lm_gg)

#########################
# V - GF                #
#########################


gf.obj <- go_gf(Y, X, X.pred, snps.set)
go.gf <- gf.obj$go
go.gf <- go.gf^2

PhenoMoyAccess <- pheno %>% 
  # Compute mean phenotypic measures per accession
  group_by(Accession) %>% 
  summarise(meanNCS=mean(NoCollSpike, na.rm=TRUE), sdNCS=sd(NoCollSpike, na.rm=TRUE),
            meanTWS=mean(TWeightSeed, na.rm=TRUE), sdTWS=sd(TWeightSeed, na.rm=TRUE),
            meanW100=mean(Weight100, na.rm=TRUE), sdW100=sd(Weight100, na.rm=TRUE),
            meanWSS=mean(TWeightSpike, na.rm=TRUE), sdWSS=sd(TWeightSpike, na.rm=TRUE)
  )
PhenoMoyAccess <- as.data.frame(PhenoMoyAccess)

df.offset.gf <- data.frame(rownames(genome_align))
colnames(df.offset.gf) <- c("Accession")
df.offset.gf["offsetgf"] <- go.gf
PhenoMoyAccess <- merge(PhenoMoyAccess, df.offset.gf)

# Plot
data_gf <- data.frame(cbind(PhenoMoyAccess$offsetgf[remove_outlier], log(PhenoMoyAccess$meanTWS[remove_outlier])))
colnames(data_gf) <- c("GF", "Fitness")

p2 <- ggplot(data_gf,aes(GF, Fitness)) +
  geom_point(shape=21, colour="black", fill="palegreen3", size=13, stroke=3) +
  geom_smooth(method='lm', se=F, col='orange1', size=6, linetype="dashed") +
  ylab("") +
  xlab("") +
  theme_bw() +
  theme(legend.position="none", 
        text = element_text(size=60))


p2

lm_gf <- lm(log(PhenoMoyAccess$meanTWS[remove_outlier]) ~ PhenoMoyAccess$offsetgf[remove_outlier])
summary(lm_gf)


#########################
# VI - GF conf          #
#########################


gf.obj <- go_gf(Y, X, X.pred, snps.set, lfmm.obj@U)
go.gf <- gf.obj$go
go.gf <- go.gf^2
PhenoMoyAccess <- pheno %>% 
  # Compute mean phenotypic measures per accession
  group_by(Accession) %>% 
  summarise(meanNCS=mean(NoCollSpike, na.rm=TRUE), sdNCS=sd(NoCollSpike, na.rm=TRUE),
            meanTWS=mean(TWeightSeed, na.rm=TRUE), sdTWS=sd(TWeightSeed, na.rm=TRUE),
            meanW100=mean(Weight100, na.rm=TRUE), sdW100=sd(Weight100, na.rm=TRUE),
            meanWSS=mean(TWeightSpike, na.rm=TRUE), sdWSS=sd(TWeightSpike, na.rm=TRUE)
  )
PhenoMoyAccess <- as.data.frame(PhenoMoyAccess)

df.offset.gf <- data.frame(rownames(genome_align))
colnames(df.offset.gf) <- c("Accession")
df.offset.gf["offsetgf"] <- go.gf
PhenoMoyAccess <- merge(PhenoMoyAccess, df.offset.gf)

# Plot
data_gf <- data.frame(cbind(PhenoMoyAccess$offsetgf[remove_outlier], log(PhenoMoyAccess$meanTWS[remove_outlier])))
colnames(data_gf) <- c("GF", "Fitness")

p2 <- ggplot(data_gf,aes(GF, Fitness)) +
  geom_point(shape=21, colour="black", fill="palegreen3", size=13, stroke=3) +
  geom_smooth(method='lm', se=F, col='orange1', size=6, linetype="dashed") +
  ylab("") +
  xlab("") +
  theme_bw() +
  theme(legend.position="none", 
        text = element_text(size=60),
        axis.text.x=element_blank(),
        axis.text.y=element_blank())


p2

lm_gf_conf <- lm(log(PhenoMoyAccess$meanTWS[remove_outlier]) ~ PhenoMoyAccess$offsetgf[remove_outlier])
summary(lm_gf)

#########################
# VI - RDA              #
#########################

rda.go <- go_rda(Y, X, X.pred, snps.set)

PhenoMoyAccess <- pheno %>% 
  # Compute mean phenotypic measures per accession
  group_by(Accession) %>% 
  summarise(meanNCS=mean(NoCollSpike, na.rm=TRUE), sdNCS=sd(NoCollSpike, na.rm=TRUE),
            meanTWS=mean(TWeightSeed, na.rm=TRUE), sdTWS=sd(TWeightSeed, na.rm=TRUE),
            meanW100=mean(Weight100, na.rm=TRUE), sdW100=sd(Weight100, na.rm=TRUE),
            meanWSS=mean(TWeightSpike, na.rm=TRUE), sdWSS=sd(TWeightSpike, na.rm=TRUE)
  )
PhenoMoyAccess <- as.data.frame(PhenoMoyAccess)

df.offset.rda <- data.frame(rownames(genome_align))
colnames(df.offset.rda) <- c("Accession")
df.offset.rda["offsetrda"] <- rda.go
PhenoMoyAccess <- merge(PhenoMoyAccess, df.offset.rda)

data_rda <- data.frame(cbind(PhenoMoyAccess$offsetrda[remove_outlier], log(PhenoMoyAccess$meanTWS[remove_outlier])))
colnames(data_rda) <- c("RDA", "Fitness")

p2 <- ggplot(data_rda,aes(RDA, Fitness)) +
  geom_point(shape=21, colour="black", fill="pink3", size=13, stroke=3) +
  geom_smooth(method='lm', se=F, col='orange1', size=6, linetype="dashed") +
  ylab("") +
  xlab("") +
  theme_bw() +
  theme(legend.position="none", 
        text = element_text(size=60),
        axis.text.x=element_blank(),
        axis.text.y=element_blank())

p2


lm_rda <- lm(log(PhenoMoyAccess$meanTWS[remove_outlier]) ~ PhenoMoyAccess$offsetrda[remove_outlier])
summary(lm_rda)


#########################
# VI - RDA CONF         #
#########################

rda.go <- go_rda_conf(Y, X, X.pred,lfmm.obj, snps.set)

PhenoMoyAccess <- pheno %>% 
  # Compute mean phenotypic measures per accession
  group_by(Accession) %>% 
  summarise(meanNCS=mean(NoCollSpike, na.rm=TRUE), sdNCS=sd(NoCollSpike, na.rm=TRUE),
            meanTWS=mean(TWeightSeed, na.rm=TRUE), sdTWS=sd(TWeightSeed, na.rm=TRUE),
            meanW100=mean(Weight100, na.rm=TRUE), sdW100=sd(Weight100, na.rm=TRUE),
            meanWSS=mean(TWeightSpike, na.rm=TRUE), sdWSS=sd(TWeightSpike, na.rm=TRUE)
  )
PhenoMoyAccess <- as.data.frame(PhenoMoyAccess)

df.offset.rda <- data.frame(rownames(genome_align))
colnames(df.offset.rda) <- c("Accession")
df.offset.rda["offsetrda"] <- rda.go
PhenoMoyAccess <- merge(PhenoMoyAccess, df.offset.rda)

data_rda <- data.frame(cbind(PhenoMoyAccess$offsetrda[remove_outlier], log(PhenoMoyAccess$meanTWS[remove_outlier])))
colnames(data_rda) <- c("RDA", "Fitness")

p2 <- ggplot(data_rda,aes(RDA, Fitness)) +
  geom_point(shape=21, colour="black", fill="pink3", size=13, stroke=3) +
  geom_smooth(method='lm', se=F, col='orange1', size=6, linetype="dashed") +
  ylab("") +
  xlab("") +
  theme_bw() +
  theme(legend.position="none", 
        text = element_text(size=60),
        axis.text.x=element_blank(),
        axis.text.y=element_blank())

p2


lm_rda_conf <- lm(log(PhenoMoyAccess$meanTWS[remove_outlier]) ~ PhenoMoyAccess$offsetrda[remove_outlier])
summary(lm_rda)


#########################
# VII - RONA            #
#########################

go.rona <- go_rona(Y, X, X.pred, snps.set)
go.rona <- go.rona^2

PhenoMoyAccess <- pheno %>% 
  # Compute mean phenotypic measures per accession
  group_by(Accession) %>% 
  summarise(meanNCS=mean(NoCollSpike, na.rm=TRUE), sdNCS=sd(NoCollSpike, na.rm=TRUE),
            meanTWS=mean(TWeightSeed, na.rm=TRUE), sdTWS=sd(TWeightSeed, na.rm=TRUE),
            meanW100=mean(Weight100, na.rm=TRUE), sdW100=sd(Weight100, na.rm=TRUE),
            meanWSS=mean(TWeightSpike, na.rm=TRUE), sdWSS=sd(TWeightSpike, na.rm=TRUE)
  )
PhenoMoyAccess <- as.data.frame(PhenoMoyAccess)

df.offset.rona <- data.frame(rownames(genome_align))
colnames(df.offset.rona) <- c("Accession")
df.offset.rona["offsetrona"] <- go.rona
PhenoMoyAccess <- merge(PhenoMoyAccess, df.offset.rona)

data_rona <- data.frame(cbind(PhenoMoyAccess$offsetrona[remove_outlier], log(PhenoMoyAccess$meanTWS[remove_outlier])))
colnames(data_rona) <- c("Rona", "Fitness")

p1 <- ggplot(data_rona,aes(Rona, Fitness)) +
  geom_point(shape=21, colour="black", fill="lightgoldenrod2", size=13, stroke=3) +
  geom_smooth(method='lm', se=F, col='orange1', size=6, linetype="dashed") +
  ylab("") +
  xlab("") +
  theme_bw() +
  theme(legend.position="none", 
        text = element_text(size=60))

p1

lm_rona <- lm(log(PhenoMoyAccess$meanTWS[remove_outlier]) ~ PhenoMoyAccess$offsetrona[remove_outlier])
summary(lm_rona)



#########################
# VIII - RONA CONf      #
#########################

go.rona <- go_rona_conf(X, X.pred, lfmm.obj, snps.set)
go.rona <- go.rona^2

PhenoMoyAccess <- pheno %>% 
  # Compute mean phenotypic measures per accession
  group_by(Accession) %>% 
  summarise(meanNCS=mean(NoCollSpike, na.rm=TRUE), sdNCS=sd(NoCollSpike, na.rm=TRUE),
            meanTWS=mean(TWeightSeed, na.rm=TRUE), sdTWS=sd(TWeightSeed, na.rm=TRUE),
            meanW100=mean(Weight100, na.rm=TRUE), sdW100=sd(Weight100, na.rm=TRUE),
            meanWSS=mean(TWeightSpike, na.rm=TRUE), sdWSS=sd(TWeightSpike, na.rm=TRUE)
  )
PhenoMoyAccess <- as.data.frame(PhenoMoyAccess)

df.offset.rona <- data.frame(rownames(genome_align))
colnames(df.offset.rona) <- c("Accession")
df.offset.rona["offsetrona"] <- go.rona
PhenoMoyAccess <- merge(PhenoMoyAccess, df.offset.rona)

data_rona <- data.frame(cbind(PhenoMoyAccess$offsetrona[remove_outlier], log(PhenoMoyAccess$meanTWS[remove_outlier])))
colnames(data_rona) <- c("Rona", "Fitness")

p1 <- ggplot(data_rona,aes(Rona, Fitness)) +
  geom_point(shape=21, colour="black", fill="lightgoldenrod2", size=13, stroke=3) +
  geom_smooth(method='lm', se=F, col='orange1', size=6, linetype="dashed") +
  ylab("") +
  xlab("") +
  xlim(0,0.13) +
  theme_bw() +
  theme(legend.position="none", 
        text = element_text(size=60),
        axis.text.x=element_blank(),
        axis.text.y=element_blank())

p1

lm_rona_conf <- lm(log(PhenoMoyAccess$meanTWS[remove_outlier]) ~ PhenoMoyAccess$offsetrona[remove_outlier])
summary(lm_rona)

#########################
# VII - Distance        #
#########################

eucl.dist <- euclidean_distance(X, X.pred)
eucl.dist <- eucl.dist^2

PhenoMoyAccess <- pheno %>% 
  # Compute mean phenotypic measures per accession
  group_by(Accession) %>% 
  summarise(meanNCS=mean(NoCollSpike, na.rm=TRUE), sdNCS=sd(NoCollSpike, na.rm=TRUE),
            meanTWS=mean(TWeightSeed, na.rm=TRUE), sdTWS=sd(TWeightSeed, na.rm=TRUE),
            meanW100=mean(Weight100, na.rm=TRUE), sdW100=sd(Weight100, na.rm=TRUE),
            meanWSS=mean(TWeightSpike, na.rm=TRUE), sdWSS=sd(TWeightSpike, na.rm=TRUE)
  )
PhenoMoyAccess <- as.data.frame(PhenoMoyAccess)

df.offset.dist <- data.frame(rownames(genome_align))
colnames(df.offset.dist) <- c("Accession")
df.offset.dist["dist"] <- eucl.dist
PhenoMoyAccess <- merge(PhenoMoyAccess, df.offset.dist)

lm_dist <- lm(log(PhenoMoyAccess$meanTWS[remove_outlier]) ~ PhenoMoyAccess$dist[remove_outlier])
summary(lm_dist)




###############################
###############################
##      All effect size      ##
###############################
###############################


snps.set <- seq(1, ncol(Y))

#########################
# IV - Genetic gap      #
#########################

genetic_gap_scale <- genetic.gap(Y, X, pred.env=X.pred, K=10, scale=T, candidate.loci=snps.set)


# relate genetic gap to log(fitness_pred)

pheno <- read.csv('./data/data_mil/RAWDATA_Pheno_ProxyYield_2016-2017_6Trials_Sadore.csv', header=T, sep=";", na.string="NA")
PhenoMoyAccess <- pheno %>% 
  # Compute mean phenotypic measures per accession
  group_by(Accession) %>% 
  summarise(meanNCS=mean(NoCollSpike, na.rm=TRUE), sdNCS=sd(NoCollSpike, na.rm=TRUE),
            meanTWS=mean(TWeightSeed, na.rm=TRUE), sdTWS=sd(TWeightSeed, na.rm=TRUE),
            meanW100=mean(Weight100, na.rm=TRUE), sdW100=sd(Weight100, na.rm=TRUE),
            meanWSS=mean(TWeightSpike, na.rm=TRUE), sdWSS=sd(TWeightSpike, na.rm=TRUE)
  )
PhenoMoyAccess <- as.data.frame(PhenoMoyAccess)
df.offset <- data.frame(rownames(genome_align))
colnames(df.offset) <- c("Accession")
df.offset["offsetcausal"] <- genetic_gap_scale$offset
PhenoMoyAccess <- merge(PhenoMoyAccess, df.offset)
remove_outlier <- log(PhenoMoyAccess$meanTWS) > 3

# And plot

data_lfmm <- data.frame(cbind(PhenoMoyAccess$offsetcausal[remove_outlier], log(PhenoMoyAccess$meanTWS[remove_outlier])))
colnames(data_lfmm) <- c("LFMM", "Fitness")


p3 <- ggplot(data_lfmm,aes(LFMM, Fitness)) +
  geom_point(shape=21, colour="black", fill="navyblue", size=13, stroke=3) +
  geom_smooth(method='lm', se=F, col='orange1', size=6, linetype="dashed") +
  ylab("") +
  xlab("") +
  theme_bw() +
  theme(legend.position="none", 
        text = element_text(size=60),
        axis.text.x=element_blank(),
        axis.text.y=element_blank())


p3

lm_gg <- lm(log(PhenoMoyAccess$meanTWS[remove_outlier]) ~ PhenoMoyAccess$offsetcausal[remove_outlier])
summary(lm_gg)


#########################
# VI - RDA              #
#########################

rda.go <- go_rda(Y, X, X.pred, snps.set)

PhenoMoyAccess <- pheno %>% 
  # Compute mean phenotypic measures per accession
  group_by(Accession) %>% 
  summarise(meanNCS=mean(NoCollSpike, na.rm=TRUE), sdNCS=sd(NoCollSpike, na.rm=TRUE),
            meanTWS=mean(TWeightSeed, na.rm=TRUE), sdTWS=sd(TWeightSeed, na.rm=TRUE),
            meanW100=mean(Weight100, na.rm=TRUE), sdW100=sd(Weight100, na.rm=TRUE),
            meanWSS=mean(TWeightSpike, na.rm=TRUE), sdWSS=sd(TWeightSpike, na.rm=TRUE)
  )
PhenoMoyAccess <- as.data.frame(PhenoMoyAccess)

df.offset.rda <- data.frame(rownames(genome_align))
colnames(df.offset.rda) <- c("Accession")
df.offset.rda["offsetrda"] <- rda.go
PhenoMoyAccess <- merge(PhenoMoyAccess, df.offset.rda)

data_rda <- data.frame(cbind(PhenoMoyAccess$offsetrda[remove_outlier], log(PhenoMoyAccess$meanTWS[remove_outlier])))
colnames(data_rda) <- c("RDA", "Fitness")

p2 <- ggplot(data_rda,aes(RDA, Fitness)) +
  geom_point(shape=21, colour="black", fill="pink3", size=13, stroke=3) +
  geom_smooth(method='lm', se=F, col='orange1', size=6, linetype="dashed") +
  ylab("") +
  xlab("") +
  theme_bw() +
  theme(legend.position="none", 
        text = element_text(size=60))

p2


lm_rda <- lm(log(PhenoMoyAccess$meanTWS[remove_outlier]) ~ PhenoMoyAccess$offsetrda[remove_outlier])
summary(lm_rda)



#########################
# VI - RDA  conf        #
#########################

rda.go <- go_rda_conf(Y, X, X.pred, lfmm.obj, snps.set)

PhenoMoyAccess <- pheno %>% 
  # Compute mean phenotypic measures per accession
  group_by(Accession) %>% 
  summarise(meanNCS=mean(NoCollSpike, na.rm=TRUE), sdNCS=sd(NoCollSpike, na.rm=TRUE),
            meanTWS=mean(TWeightSeed, na.rm=TRUE), sdTWS=sd(TWeightSeed, na.rm=TRUE),
            meanW100=mean(Weight100, na.rm=TRUE), sdW100=sd(Weight100, na.rm=TRUE),
            meanWSS=mean(TWeightSpike, na.rm=TRUE), sdWSS=sd(TWeightSpike, na.rm=TRUE)
  )
PhenoMoyAccess <- as.data.frame(PhenoMoyAccess)

df.offset.rda <- data.frame(rownames(genome_align))
colnames(df.offset.rda) <- c("Accession")
df.offset.rda["offsetrda"] <- rda.go
PhenoMoyAccess <- merge(PhenoMoyAccess, df.offset.rda)

data_rda <- data.frame(cbind(PhenoMoyAccess$offsetrda[remove_outlier], log(PhenoMoyAccess$meanTWS[remove_outlier])))
colnames(data_rda) <- c("RDA", "Fitness")

p2 <- ggplot(data_rda,aes(RDA, Fitness)) +
  geom_point(shape=21, colour="black", fill="pink4", size=13, stroke=3) +
  geom_smooth(method='lm', se=F, col='orange1', size=6, linetype="dashed") +
  ylab("") +
  xlab("") +
  theme_bw() +
  theme(legend.position="none", 
        text = element_text(size=60),
        axis.text.x=element_blank(),
        axis.text.y=element_blank())

p2


lm_rda <- lm(log(PhenoMoyAccess$meanTWS[remove_outlier]) ~ PhenoMoyAccess$offsetrda[remove_outlier])
summary(lm_rda)


#########################
# VII - RONA            #
#########################

go.rona <- go_rona(Y, X, X.pred, snps.set)
go.rona <- go.rona^2

PhenoMoyAccess <- pheno %>% 
  # Compute mean phenotypic measures per accession
  group_by(Accession) %>% 
  summarise(meanNCS=mean(NoCollSpike, na.rm=TRUE), sdNCS=sd(NoCollSpike, na.rm=TRUE),
            meanTWS=mean(TWeightSeed, na.rm=TRUE), sdTWS=sd(TWeightSeed, na.rm=TRUE),
            meanW100=mean(Weight100, na.rm=TRUE), sdW100=sd(Weight100, na.rm=TRUE),
            meanWSS=mean(TWeightSpike, na.rm=TRUE), sdWSS=sd(TWeightSpike, na.rm=TRUE)
  )
PhenoMoyAccess <- as.data.frame(PhenoMoyAccess)

df.offset.rona <- data.frame(rownames(genome_align))
colnames(df.offset.rona) <- c("Accession")
df.offset.rona["offsetrona"] <- go.rona
PhenoMoyAccess <- merge(PhenoMoyAccess, df.offset.rona)

data_rona <- data.frame(cbind(PhenoMoyAccess$offsetrona[remove_outlier], log(PhenoMoyAccess$meanTWS[remove_outlier])))
colnames(data_rona) <- c("Rona", "Fitness")

p1 <- ggplot(data_rona,aes(Rona, Fitness)) +
  geom_point(shape=21, colour="black", fill="lightgoldenrod2", size=13, stroke=3) +
  geom_smooth(method='lm', se=F, col='orange1', size=6, linetype="dashed") +
  ylab("") +
  xlab("") +
  theme_bw() +
  theme(legend.position="none", 
        text = element_text(size=60))

p1

lm_rona <- lm(log(PhenoMoyAccess$meanTWS[remove_outlier]) ~ PhenoMoyAccess$offsetrona[remove_outlier])
summary(lm_rona)

#########################
# VII - RONA  conf      #
#########################

go.rona <- go_rona_conf( X, X.pred, lfmm.obj, snps.set)
go.rona <- go.rona^2

PhenoMoyAccess <- pheno %>% 
  # Compute mean phenotypic measures per accession
  group_by(Accession) %>% 
  summarise(meanNCS=mean(NoCollSpike, na.rm=TRUE), sdNCS=sd(NoCollSpike, na.rm=TRUE),
            meanTWS=mean(TWeightSeed, na.rm=TRUE), sdTWS=sd(TWeightSeed, na.rm=TRUE),
            meanW100=mean(Weight100, na.rm=TRUE), sdW100=sd(Weight100, na.rm=TRUE),
            meanWSS=mean(TWeightSpike, na.rm=TRUE), sdWSS=sd(TWeightSpike, na.rm=TRUE)
  )
PhenoMoyAccess <- as.data.frame(PhenoMoyAccess)

df.offset.rona <- data.frame(rownames(genome_align))
colnames(df.offset.rona) <- c("Accession")
df.offset.rona["offsetrona"] <- go.rona
PhenoMoyAccess <- merge(PhenoMoyAccess, df.offset.rona)

data_rona <- data.frame(cbind(PhenoMoyAccess$offsetrona[remove_outlier], log(PhenoMoyAccess$meanTWS[remove_outlier])))
colnames(data_rona) <- c("Rona", "Fitness")

p1 <- ggplot(data_rona,aes(Rona, Fitness)) +
  geom_point(shape=21, colour="black", fill="lightgoldenrod2", size=13, stroke=3) +
  geom_smooth(method='lm', se=F, col='orange1', size=6, linetype="dashed") +
  ylab("") +
  xlab("") +
  theme_bw() +
  theme(legend.position="none", 
        text = element_text(size=60))

p1

lm_rona <- lm(log(PhenoMoyAccess$meanTWS[remove_outlier]) ~ PhenoMoyAccess$offsetrona[remove_outlier])
summary(lm_rona)

###############################
###############################
##     Combining offset     ##
###############################
###############################


