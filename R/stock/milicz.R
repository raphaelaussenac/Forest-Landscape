###############################################################
# initialisation
###############################################################

# clean up environment
rm(list = ls())

# print options
options(digits=10)

# load packages
library(tidyverse)

# set work directory
setwd('C:/Users/raphael.aussenac/Documents/GitHub/LandscapeInit')

# function to replace commas by dots and transform character vectors
# into numeric vectors
comTodot <- function(charVect){
  charVect <- as.numeric(str_replace_all(charVect, ',', '\\.'))
  return(charVect)
}

# load data

# plot data
plotdf <- read.csv('./data/milicz/Milicz_field_plots.csv', sep = ';')
plotdf$X <- comTodot(plotdf$X)
plotdf$Y <- comTodot(plotdf$Y)
plotdf$BA <- comTodot(plotdf$BA)
plotdf$QM_DBH <- comTodot(plotdf$QM_DBH)
plotdf$decid_prop <- comTodot(plotdf$decid_prop)

# tree data
treedf <- read.csv('./data/milicz/Milicz_trees.csv', sep = ';')
treedf$DBH <- comTodot(treedf$DBH)
treedf$height <- comTodot(treedf$height)
treedf$species <- as.factor(treedf$species)
treedf$species <- recode_factor(treedf$species, 'Malus silvestris' = 'Malus sylvestris',
                                                'Quercus undefined' = 'Quercus spp.',
                                                'Rhamnus frangula' = 'Frangula alnus',
                                                'Robinia pseudoacacia' = 'Robinia pseudacacia',
                                                'Ulmus' = 'Ulmus spp.')
#
