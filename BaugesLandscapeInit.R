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
names(elevation) <- 'elev'

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
names(parkRaster) <- 'park'

# convert forest into raster and set extent + resolution
forestRaster <- rasterize(forest, r, field = "ID")
# set projection
crs(forestRaster) <- crs(forest)
# convert NA into 0
isBecomes <- cbind(c(1:nrow(forest), NA),
                   c(rep(1, nrow(forest)), 0))
forestRaster <- reclassify(forestRaster, rcl = isBecomes)
names(forestRaster) <- 'forest'

# swhc (cm)
swhc <- raster("./data/GEO/rum_500_v2009.tif")
swhc <- resample(swhc, elevation)
swhc <- swhc/10 # convert into cm
names(swhc) <- 'swhc'

# Quadratic diameter (cm) [0, 80]
dg <- raster('./data/GEO/rastDg75error.clean.tif')
crs(dg) <- crs(park)
dg <- resample(dg, elevation)

# basal area (m2) [0, 120]
BA <- raster('./data/GEO/rastG75error.clean.tif')
crs(BA) <- crs(park)
BA <- resample(BA, elevation)

# # Quadratic diameter (cm) with all extreme values
# dgExtrem <- raster('./data/GEO/rastDg75error.forest.tif')
# crs(dgExtrem) <- crs(park)
# dgExtrem <- resample(dgExtrem, elevation)
#
# # basal area (m2) with all extreme values
# BAExtrem <- raster('./data/GEO/rastG75error.forest.tif')
# crs(BAExtrem) <- crs(park)
# BAExtrem <- resample(BAExtrem, elevation)

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

# create cell ID raster
ext <- floor(extent(elevation))
cellID <- raster(ext, res=res(elevation))
cellID$cellID <- c(1:(nrow(cellID) * ncol(cellID)))
cellID <- cellID$cellID

# greco
greco <- readOGR(dsn = "./data/GEO", layer = "greco_l93", encoding = "UTF-8", use_iconv = TRUE)
greco <- spTransform(greco, crs(park)) # change projection
# convert into a raster
grecoRaster <- rasterize(greco, cellID, field="CODEGRECO")
crs(grecoRaster) <- crs(park)
names(grecoRaster) <- 'GRECO'

# geol
geol <- readOGR(dsn = "./data/GEO", layer = "geol", encoding = "UTF-8", use_iconv = TRUE)
geol <- geol[, 'NOTATION']
classGeol <- read.csv("./data/geol/classificationGeol.csv", header = TRUE, sep = ";")
classGeol <- classGeol[, c('NOTATION', 'Code_carbonate', 'Code_hydro')]
colnames(classGeol) <- c('NOTATION', 'Cd_crbn', 'Cd_hydr')
geol <- merge(geol, classGeol, by = 'NOTATION')
# convert into a raster
Cd_crbn <- rasterize(geol, cellID, field = "Cd_crbn")
crs(Cd_crbn) <- crs(park)
names(Cd_crbn) <- 'Cd_crbn'
Cd_hydr <- rasterize(geol, cellID, field = "Cd_hydr")
crs(Cd_hydr) <- crs(park)
names(Cd_hydr) <- 'Cd_hydr'

###############################################################
# save ascii
###############################################################

# create directory to save plots
if (!(dir.exists('./data/init'))) {dir.create('./data/init', recursive = TRUE)}
if (!(dir.exists('./initialLandscape'))) {dir.create('./initialLandscape', recursive = TRUE)}

writeRaster(parkRaster, filename = "./initialLandscape/parkMask.asc", format = "ascii", overwrite = TRUE)
writeRaster(forestRaster, filename = "./initialLandscape/forestMask.asc", format = "ascii", overwrite = TRUE)
writeRaster(elevation, filename = "./initialLandscape/elev.asc", format = "ascii", overwrite = TRUE)
writeRaster(slope, filename = "./initialLandscape/slope.asc", format = "ascii", overwrite = TRUE)
writeRaster(aspect, filename = "./initialLandscape/aspect.asc", format = "ascii", overwrite = TRUE)
writeRaster(swhc, filename = "./initialLandscape/swhc.asc", format = "ascii", overwrite = TRUE)
writeRaster(dg, filename = "./data/init/dg.asc", format = "ascii", overwrite = TRUE)
writeRaster(BA, filename = "./data/init/BA.asc", format = "ascii", overwrite = TRUE)
writeRaster(N, filename = "./data/init/N.asc", format = "ascii", overwrite = TRUE)
writeRaster(Dprop, filename = "./data/init/Dprop.asc", format = "ascii", overwrite = TRUE)
writeRaster(cellID, filename = "./data/init/cellID.asc", format = "ascii", overwrite = TRUE)
writeRaster(grecoRaster, filename = "./data/init/greco.asc", format = "ascii", overwrite = TRUE)
writeRaster(Cd_crbn, filename = "./data/init/Cd_crbn.asc", format = "ascii", overwrite = TRUE)
writeRaster(Cd_hydr, filename = "./data/init/Cd_hydr.asc", format = "ascii", overwrite = TRUE)

###############################################################
# save data frame
###############################################################

# create raster stack
rasterStack <- stack(cellID, parkRaster, forestRaster, elevation, slope,
                     aspect, swhc, grecoRaster, Cd_crbn, Cd_hydr)
# plot(rasterStack)

# convert into data frame
envdf <- as.data.frame(rasterStack)

# rename back GRECO regions
envdf[envdf$GRECO == 4, 'GRECO'] <- 'C'
envdf[envdf$GRECO == 6, 'GRECO'] <- 'E'
envdf[envdf$GRECO == 9, 'GRECO'] <- 'H'

# save
write.csv(envdf, file = './initialLandscape/envVariables.csv', row.names = FALSE)
