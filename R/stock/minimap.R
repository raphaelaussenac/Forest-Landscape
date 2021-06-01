###############################################################
# initialisation
###############################################################

# clean up environment
rm(list = ls())

# load packages
library(raster)
library(dplyr)

# set work directory
setwd('C:/Users/raphael.aussenac/Documents/GitHub/LandscapeInit')

###############################################################
# 1 - open all rasters
###############################################################

aspect <- raster('./initialLandscape/aspect.asc')
cellID25 <- raster('./initialLandscape/cellID25.asc')
cellID100 <- raster('./initialLandscape/cellID100.asc')
elev <- raster('./initialLandscape/elev.asc')
forest <- raster('./initialLandscape/forestMask.asc')
park <- raster('./initialLandscape/parkMask.asc')
pH <- raster('./initialLandscape/pH.asc')
slope <- raster('./initialLandscape/slope.asc')
swhc <- raster('./initialLandscape/swhc.asc')

###############################################################
# 2 - crop on smaller extent
###############################################################

ext <- extent(945000, 945500, 6513400, 6513900)
aspect <- crop(aspect, ext)
cellID25 <- crop(cellID25, ext)
cellID100 <- crop(cellID100, ext, snap = 'out')
elev <- crop(elev, ext)
forest <- crop(forest, ext)
park <- crop(park, ext)
pH <- crop(pH, ext)
slope <- crop(slope, ext)
swhc <- crop(swhc, ext)

###############################################################
# 3 - save rasters
###############################################################

writeRaster(park, filename = './initialLandscape/minimap/parkMask.asc', format = 'ascii', overwrite = TRUE)
writeRaster(forest, filename = './initialLandscape/minimap/forestMask.asc', format = 'ascii', overwrite = TRUE)
writeRaster(elev, filename = './initialLandscape/minimap/elev.asc', format = 'ascii', overwrite = TRUE)
writeRaster(slope, filename = './initialLandscape/minimap/slope.asc', format = 'ascii', overwrite = TRUE)
writeRaster(aspect, filename = './initialLandscape/minimap/aspect.asc', format = 'ascii', overwrite = TRUE)
writeRaster(swhc, filename = './initialLandscape/minimap/swhc.asc', format = 'ascii', overwrite = TRUE)
writeRaster(pH, filename = './initialLandscape/minimap/pH.asc', format = 'ascii', overwrite = TRUE)
writeRaster(cellID25, filename = './initialLandscape/minimap/cellID25.asc', format = 'ascii', overwrite = TRUE)
writeRaster(cellID100, filename = './initialLandscape/minimap/cellID100.asc', format = 'ascii', overwrite = TRUE)

###############################################################
# 3 - transform tree and environmental data frames
###############################################################

# retrieve cellID
cellID <- values(cellID25)

# reduce envVariables df
df <- read.csv('./initialLandscape/envVariables.csv')
df <- df[df$cellID25 %in% cellID,]
# calculate nb of forest cells per ha for the mini landscape
# may be different from nb of forest cells in envVariables.csv
df <- df %>% group_by(cellID100) %>% mutate(forestCellsPerHa = sum(forest))
# plot hist of nb of forest cells per ha
test <- df %>% group_by(cellID100) %>% summarise(forestCellsPerHa = unique(forestCellsPerHa))
hist(test$forestCellsPerHa, breaks = c(0:16))
# save
write.csv(df, './initialLandscape/minimap/envVariables.csv', row.names = F)

# surface (ha)
nrow(df)*625/10000

# reduce trees75 df
df <- read.csv('./initialLandscape/trees75.csv')
df <- df[df$cellID25 %in% cellID,]
write.csv(df, './initialLandscape/minimap/trees75.csv', row.names = F)
