###############################################################
# initialisation
###############################################################

# clean up environment
rm(list = ls())

# load packages
library(rgdal)
library(raster)
library(velox)
library(gdalUtils)
library(sf)

# set work directory
setwd("C:/Users/raphael.aussenac/Documents/GitHub/LandscapeInit")

###############################################################
# load/create data and set extent and resolution
###############################################################

# park
park <- readOGR(dsn = "./data/GEO", layer = "park", encoding = "UTF-8", use_iconv = TRUE)

# forest
forest <- readOGR(dsn = "./data/GEO", layer = "BD_Foret_V2_PNRfilled_Foret_2014", encoding = "UTF-8", use_iconv = TRUE)

# elevation (m a.s.l.)
elevation <- raster("./data/GEO/MNT_all_5m.tif")
# set projection
crs(elevation) <- crs(park)
# set resolution
elevation <- aggregate(elevation, fact = 5)

# slope (degree)
slope <- terrain(elevation, opt = 'slope', unit = 'degrees', neighbors = 8)

# aspect (degree)
aspect <- terrain(elevation, opt = 'aspect', unit = 'degrees', neighbors = 8)

# convert park into raster and set extent + resolution
ext <- floor(extent(elevation))
r <- raster(ext, res=res(elevation))
parkRaster <- rasterize(park, r, field = "ID")
# set projection
crs(parkRaster) <- crs(park)
# convert NA into 0
isBecomes <- cbind(c(1, NA),
                   c(1, 0))
parkRaster <- reclassify(parkRaster, rcl = isBecomes)

# convert forest into raster and set extent + resolution
forestRaster <- rasterize(forest, r, field = "ID")
# set projection
crs(forestRaster) <- crs(forest)
# convert NA into 0
isBecomes <- cbind(c(1:nrow(forest), NA),
                   c(rep(1, nrow(forest)), 0))
forestRaster <- reclassify(forestRaster, rcl = isBecomes)

# swhc (cm)
swhc <- raster("./data/GEO/rum_500_v2009.tif")
swhc <- resample(swhc, elevation)
swhc <- swhc/10 # convert into cm

# Quadratic diameter (cm)
dg <- raster('./data/GEO/rastDg75error.clean.tif')
crs(dg) <- crs(park)
dg <- resample(dg, elevation)

# basal area (m2)
BA <- raster('./data/GEO/rastG75error.clean.tif')
crs(BA) <- crs(park)
BA <- resample(BA, elevation)

# N
N <- raster('./data/GEO/N_pred.tif')
crs(N) <- crs(park)
N <- resample(N, elevation)

# large tree (BA or proportion of total BA?)
# LTBA <- raster('./data/GEO/GGB_pred.tif')
# crs(LTBA) <- crs(park)
# LTBA <- resample(LTBA, elevation)

# Deciduous proportion (% of total BA)
Dprop <- raster('./data/GEO/propGR_ONF_25filled.tif')
Dprop <- 100 - Dprop
crs(Dprop) <- crs(park)
Dprop <- resample(Dprop, elevation)

###############################################################
# save ascii
###############################################################

# create directory to save plots
if (!(dir.exists('Init'))) {dir.create('Init', recursive = TRUE)}

writeRaster(parkRaster, filename = "./Init/parkMask.asc", format = "ascii", datatype = 'INT4S', overwrite = TRUE)
writeRaster(forestRaster, filename = "./Init/forestMask.asc", format = "ascii", datatype = 'INT4S', overwrite = TRUE)
writeRaster(elevation, filename = "./Init/elev.asc", format = "ascii", datatype = 'INT4S', overwrite = TRUE)
writeRaster(slope, filename = "./Init/slope.asc", format = "ascii", datatype = 'INT4S', overwrite = TRUE)
writeRaster(aspect, filename = "./Init/aspect.asc", format = "ascii", datatype = 'INT4S', overwrite = TRUE)
writeRaster(swhc, filename = "./Init/swhc.asc", format = "ascii", datatype = 'INT4S', overwrite = TRUE)
writeRaster(dg, filename = "./Init/dg.asc", format = "ascii", datatype = 'INT4S', overwrite = TRUE)
writeRaster(BA, filename = "./Init/BA.asc", format = "ascii", datatype = 'INT4S', overwrite = TRUE)
writeRaster(N, filename = "./Init/N.asc", format = "ascii", datatype = 'INT4S', overwrite = TRUE)
# writeRaster(LTBA, filename = "./Init/LTBA.asc", format = "ascii", datatype = 'INT4S', overwrite = TRUE)
writeRaster(Dprop, filename = "./Init/Dprop.asc", format = "ascii", datatype = 'INT4S', overwrite = TRUE)
