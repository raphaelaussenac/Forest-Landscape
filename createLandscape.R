###############################################################
# initialisation
###############################################################
# Rprof()
# clean up environment
rm(list = ls())

# set work directory
setwd('~/Documents/code/Forest-Landscape')

# load source (all R files in R folder)
file.sources = list.files('./R', pattern = '*.R', full.names = TRUE)
sapply(file.sources, source, .GlobalEnv)

###############################################################
# select landscape (bauges, milicz, sneznik)
# set number of cores
###############################################################

landscape <- 'bauges'
cores <- 6

###############################################################
# create virtual landscape
###############################################################

# running time
start_time_gen <- Sys.time()

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
start_time_compo <- Sys.time()
compoNew(landscape, cores)
end_time_compo <- Sys.time()
print(end_time_compo - start_time_compo)

# evaluate composition
# evalCompo()

# assign individual trees to all 25*25m cells
start_time_dendro <- Sys.time()
dendroNew(landscape, cores)
end_time_dendro <- Sys.time()
print(end_time_dendro - start_time_dendro)
# summaryRprof()$by.total
# summaryRprof()$by.self
# Rprof(NULL)

# produce landscape tree and environmental data
start_time_save <- Sys.time()
saveLandscape(landscape)
end_time_save <- Sys.time()
print(end_time_save - start_time_save)

# evaluate dendro variables
# evalDendro()

# model tree height
start_time_height <- Sys.time()
heightPred(landscape)
end_time_height <- Sys.time()
print(end_time_height - start_time_height)

# compare tree heights against LiDAR heights
# evalHeight(landscape)

# # create all alternative management table
# altManagTable(landscape)

# # synthesise management of all alternative scenarios
# altManagSynth(landscape)

# # evaluate scenarios
# evalSce(landscape)

# running time
end_time_gen <- Sys.time()
print(end_time_gen - start_time_gen)
