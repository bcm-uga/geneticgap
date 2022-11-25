# A quantitative theory for genomic offset statistics

This repository stores all the relevant R script to generate figure for our papers : "A quantitative theory for genomic offset statistics". The required packages are :
1. LEA
2. qvalue
3. gradientForest
4. ggplot2
5. cowplot
6. fields
7. lmtest
8. dplyr
9. tidyverse
10. rnaturalearth
11. geosphere
12. RColorBrewer



In order to run the scripts, you'll need to download the data and to put the **data** directory at the root of the repository.

## Organisation of the repository

### Data

The data repository contains 4 folders : 
1. **expfit_2variables** : It corresponds to the data generate with SLiM for the scenario "high confounding weakly polygenic" (hcwp)
2. **expfit_2variables_poly** : It corresponds to the data generate with SLiM for the scenario "high confounding highly polygenic" (hchp)
3. **poly_exp** : It corresponds to the data generate with SLiM for the scenario "low confounding highly polygenic" (lchp)
4. **poly_small** : It corresponds to the data generate with SLiM for the scenario "low confounding weakly polygenic" (lchp)
5. **data_mil** : All the data related to the pearl millet experiment

### R

At the root of the R repository, we find 2 R files :
1. **offset.R** : It contains all the function that are useful to compute offset and environmental distance
2. **utils.R** : It contains relevant function mainly to extract and transform the data from Data repository

It also contains 4 folders. In each folder, we can find scripts to obtain the principal figures and all the related supplementary figures

### Results

In order not to have to compute all the offset each time we want to plot the figures, we stored all the results in files and we stored these results in the **results** repository

### slimwork

All the slim scripts that generated the simulated data