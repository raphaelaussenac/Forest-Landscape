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

###############################################################
# select landscape (bauges, milicz)
###############################################################

landscape <- 'bauges'

###############################################################
# create virtual landscape
###############################################################

# define path to temp files
tempPath <- paste0('./data/temp/', landscape)
landPath <- paste0('./', landscape, 'Landscape')
evalPath <- paste0(landPath, '/evaluation')

# create temp, landscape and evaluation folders
if (!(dir.exists(tempPath))) {dir.create(tempPath, recursive = TRUE)}
if (!(dir.exists(landPath))) {dir.create(landPath, recursive = TRUE)}
if (!(dir.exists(evalPath))) {dir.create(evalPath, recursive = TRUE)}

# Prepare landscape data
if(landscape == 'bauges'){
  prepBauges()
}

# assign composition to all 25*25m cells
compo(landscape)

# evaluate composition
evalCompo()

# assign individual trees to all 25*25m cells
# dendro(landscape)
