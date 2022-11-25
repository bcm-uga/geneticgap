library(dplyr)



setwd('~/Documents/Rwork/diversityoffset/')


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

climate_temperature <- climate_only_mil[, 67:150]
pc_temp <- prcomp(climate_temperature)
X.temp <- pc_temp$x

genome_align <- as.matrix(t(genome_align))


pc_mil <- prcomp(climate_only_mil)
X <- pc_mil$x

PixelID.Sadore <- "39_15"
Pixels.value <- read.table("./data/data_mil/ClimateMetrics_pixelsINarea_RealObservation.txt",
                           header=TRUE, sep="\t", 
                           fileEncoding="iso-8859-1")

sadore.env <- Pixels.value[Pixels.value$PixelID==PixelID.Sadore,][4:177]
sadore.env <- sadore.env[, colsd != 0]
sadore.temp <- sadore.env[, 67:150]
X.temp.sadore <- predict(pc_temp, sadore.temp)


X.sadore <- predict(pc_mil, sadore.env)

dif.env <- X[,1] - X.sadore[1]
dif.temp <- X.temp[,1] - X.temp.sadore[1]

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
df.offset[,"env_dist"] <- dif.env
df.offset[,"temp_dist"] <- dif.temp

PhenoMoyAccess <- merge(PhenoMoyAccess, df.offset)

remove_outlier <- log(PhenoMoyAccess$meanTWS) > 3

df_s11 <- data.frame(cbind(PhenoMoyAccess$env_dist[remove_outlier], log(PhenoMoyAccess$meanTWS[remove_outlier])))
colnames(df_s11) <- c("Dif", "Fitness")

ggplot(df_s11,aes(Dif, Fitness)) +
  geom_point(shape=21, colour="black", fill="turquoise", size=13, stroke=3) +
  geom_smooth(method='loess', se=F, col='tomato', size=6) +
  ylab("") +
  xlab("") +
  theme_bw() +
  theme(legend.position="none", 
        text = element_text(size=60),
        axis.text.x=element_blank(),
        axis.text.y=element_blank())
