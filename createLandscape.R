###############################################################
# initialisation
###############################################################

# clean up environment
rm(list = ls())

# set work directory
setwd('./Documents/code/Forest-Landscape')

# load source (all R files in R folder)
file.sources = list.files('./R', pattern = '*.R', full.names = TRUE)
sapply(file.sources, source, .GlobalEnv)

###############################################################
# select landscape (bauges, milicz)
###############################################################

landscape <- 'bauges'

###############################################################
# create virtual landscape
###############################################################

# running time
start_time <- Sys.time()

# define folder structure
tempPath <- paste0('./temp/', landscape)
landPath <- paste0('./', landscape)
evalPath <- paste0(landPath, '/evaluation')
evalHeightPath <- './evalHeight'

# create folders
if (!(dir.exists(tempPath))) {dir.create(tempPath, recursive = TRUE)}
if (!(dir.exists(landPath))) {dir.create(landPath, recursive = TRUE)}
if (!(dir.exists(evalPath))) {dir.create(evalPath, recursive = TRUE)}
if (!(dir.exists(evalHeightPath))) {dir.create(evalHeightPath, recursive = TRUE)}

# Prepare landscape data
if(landscape == 'bauges'){
  prepBauges()
} else if (landscape == 'milicz'){
  prepMilicz()
} else if (landscape == 'sneznik'){
  prepSneznik()
}

# assign composition to all 25*25m cells
compo(landscape)

# evaluate composition
evalCompo()

# assign individual trees to all 25*25m cells
dendro()

# produce landscape tree and environmental data
saveLandscape(landscape)

# evaluate dendro variables
evalDendro()

# model tree height
heightPred(landscape)

# compare tree heights against LiDAR heights
evalHeight(landscape)

# create all alternative management table
altManagTable(landscape)

# synthesise management of all alternative scenarios
altManagSynth(landscape)

# evaluate scenarios
evalSce(landscape)

# running time
end_time <- Sys.time()
end_time - start_time
