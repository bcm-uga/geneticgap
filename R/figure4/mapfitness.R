library(ggplot2)
library(qvalue)
library(LEA)
library(dplyr)
library(tidyverse)
library(rnaturalearth)
library(geosphere)

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

snps.set <- seq(1,ncol(Y))
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

mergedata <- merge(PhenoMoyAccess, climate_mil)

latlongfitness <- mergedata[,c("Longitude.DD", "Latitude.DD","meanTWS")]

world <- map_data("world")



# We need to create a raster of data.
# We only have a few points so we have to create a regular grid 

Pixels.historicalMetrics <- 
  read.table("./data/data_mil/ClimateMetrics_pixelsINarea_RealObservation.txt",
             header=TRUE, sep="\t", 
             fileEncoding="iso-8859-1")

PixelID.Sadore <- "39_15"
longlat_sadore <- Pixels.historicalMetrics[Pixels.historicalMetrics$PixelID==PixelID.Sadore, c("Longitude.DD", "Latitude.DD")]
df_sadore <- as.data.frame(longlat_sadore)
colnames(df_sadore) <- c("Long", "Lat")


longlat <- Pixels.historicalMetrics[,c("Longitude.DD", "Latitude.DD")]


# We need to add points

order_longitude <- unique(sort(longlat$Longitude.DD))
new_points <- c()
for (long in order_longitude){
  min_lat <- min(longlat[longlat$Longitude.DD==long, "Latitude.DD"])
  max_lat <- max(longlat[longlat$Longitude.DD==long, "Latitude.DD"])
  current_points <- rbind(c(long, min_lat - 0.5), c(long, min_lat - 1), c(long, max_lat + 0.5), c(long, max_lat + 1))
  new_points <- rbind(new_points, current_points)
}

order_latitude <- unique(sort(longlat$Latitude.DD))
for (lat in order_latitude){
  max_long <- max(longlat[longlat$Latitude.DD==lat, "Longitude.DD"])
  current_points <- rbind(c(max_long+0.5,lat ), c(max_long + 1, lat))
  new_points <- rbind(new_points, current_points)
}

colnames(new_points) <- colnames(longlat)
longlat <- rbind(longlat, new_points)

nb_longlat <- nrow(longlat)
nb_fitness <- nrow(latlongfitness)



######

list_fitness <- c()

for (i in seq(1, nb_longlat)){
  longlat1 <- longlat[i, c("Longitude.DD", "Latitude.DD")]
  list_dist <- c()
  for (j in seq(1, nb_fitness)){
    list_dist <- c(list_dist, distm(longlat1, latlongfitness[j , c("Longitude.DD", "Latitude.DD")], fun=distHaversine))
  }
  
  index_min <- which.min(list_dist)
  list_fitness <- c(list_fitness, latlongfitness[index_min, 3])
}

df_fitness <- cbind(longlat, list_fitness)
df_fitness <- as.data.frame(df_fitness)
colnames(df_fitness) <- c("Long", "Lat", "meanTWS")


ggplot() +
  geom_map(
    data = world, map = world,
    aes(long, lat, map_id = region),
    color = "black", fill = "white", size = 0.1
  ) +
  xlim(-20,30)+
  ylim(0,20) +
  geom_raster(data = df_fitness, aes(x = Long, y = Lat, fill = meanTWS)) +
  geom_point(
    data = latlongfitness,
    aes(Longitude.DD, Latitude.DD), color = "black"
  ) +
  scale_fill_gradient(low="blue", high="red", name='') +
  geom_point(
    data = df_sadore,
    aes(Long, Lat), shape=21, color = "white", fill="black", size=4, stroke=2
  ) 

