###############################################################
# initialisation
###############################################################

# clean up environment
rm(list = ls())

# set work directory
setwd('C:/Users/raphael.aussenac/Documents/GitHub/LandscapeInit')

# load source (all R files in R folder)
# file.sources = list.files('./R', pattern = '*.R', full.names = TRUE)
# sapply(file.sources, source, .GlobalEnv)
source('./R/baugesLandscapePrep.R')
source('./R/salemSIpred.R')
source('./R/compo.R')
source('./R/spTransform.R')
source('./R/evalCompo.R')
source('./R/dendro.R')
source('./R/saveLandscape.R')

source('./R/minimap.R')

###############################################################
# select landscape (bauges, milicz)
###############################################################

landscape <- 'bauges'

###############################################################
# create virtual landscape
###############################################################

# define folder structure
tempPath <- paste0('./data/temp/', landscape)
landPath <- paste0('./', landscape, 'Landscape')
evalPath <- paste0(landPath, '/evaluation')
miniPath <- paste0(landPath, '/minimap')

# create temp, landscape and evaluation folders
if (!(dir.exists(tempPath))) {dir.create(tempPath, recursive = TRUE)}
if (!(dir.exists(landPath))) {dir.create(landPath, recursive = TRUE)}
if (!(dir.exists(evalPath))) {dir.create(evalPath, recursive = TRUE)}
if (!(dir.exists(miniPath))) {dir.create(miniPath, recursive = TRUE)}

# Prepare landscape data
if(landscape == 'bauges'){
  prepBauges()
}

# assign composition to all 25*25m cells
compo(landscape)

# evaluate composition
evalCompo()

# assign individual trees to all 25*25m cells
dendro(landscape)

# produce landscape tree and environmental data
saveLandscape()

# evaluate dendro variables
# evalDendro()

# produce minimap
minimap(landscape)
