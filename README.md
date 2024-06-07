# A quantitative theory for genomic offset statistics

This repository stores all relevant R scripts to generate figures for the paper : "A quantitative theory for genomic offset statistics". It also contains a tutorial on  how to display predictions of genomic offset statistics in geographic maps. 

## Tutorial

The tutorial (Rmd and html) is a quite general example that extracts data from the worldclim database and from IPPC sixth report future climate scenarios. It shows how to display genomic offset predictions in a geographic map using R packages for spatial analysis.  The tutorial is not related to the Pearl millet data analysis done in our study, and for which the data were not based on worldclim variables. 


## Organisation of the repository

### Data

The data repository contains 4 folders: 
1. **expfit_2variables** corresponds to the data generated with SLiM for the scenario "high confounding weakly polygenic" (hcwp)
2. **expfit_2variables_poly** corresponds to the data generated with SLiM for the scenario "high confounding highly polygenic" (hchp)
3. **poly_exp** corresponds to the data generated with SLiM for the scenario "low confounding highly polygenic" (lchp)
4. **poly_small** corresponds to the data generated with SLiM for the scenario "low confounding weakly polygenic" (lchp)
5. **data_mil** corresponds to data related to the pearl millet experiment (see Rhone et al. 2020 for access to the source data set)

In order to run the scripts,  the **data** directory must be at the root of the repository.


### R

At the root of the R repository,
1. **offset.R** contains all the functions that are useful to compute offset and environmental distances
2. **utils.R** contains functions to extract and transform data from Data repository

The repository contains 4 folders. Each folder contains scripts related to the main document figures and related supplementary figures

### Results

All results were stored in separated files

### slimwork

All SLiM scripts that generated the simulated data


## R packages required

For genotype-environment association studies and genetic offset calculations,
1. LEA 
2. gradientForest

For spatial analysis,
3. terra
4. fields
5. geodata
6. maps
7. rnaturalearth
8. geosphere

For statistics,
9. lmtest
10. qvalue

For general data analysis,
11. ggplot2
12. cowplot
13. dplyr
14. tidyverse
15. RColorBrewer

