prepSneznik <- function(){

  ###############################################################
  # initialisation
  ###############################################################

  # load packages
  require(rgdal)
  require(raster)
  require(sf)
  require(stringr)
  require(dplyr)
  require(readxl)

  ###############################################################
  # load tree data
  ###############################################################

  # load tree inventory data
  tree <- read_excel('./data/sneznik/inventory/Sneznik_tree_data.xlsx', sheet = 1)

  # replace plot id with letters by numeric ids
  tree$PLOTID <- as.numeric(as.factor(tree$PLOTID))

  # rename undefined species
  tree$SPECIES <- as.factor(tree$SPECIES)
  tree$SPECIES <- recode_factor(tree$SPECIES, 'Tilia sp.' = 'Tilia platyphyllos',
                                              'Salix sp.' = 'Salix caprea',
                                              'Pinus sp.' = 'Pinus sylvestris')

  # add missing columns remove useless columns
  tree <- tree %>% dplyr::select(PLOTID, DBH, SPECIES) %>%
                   rename(idp = PLOTID, species_name = SPECIES) %>%
                   mutate(idp = as.factor(idp), w = NA, spType = NA, CODE_TFV = as.factor(99)) %>%
                   relocate(idp, w, CODE_TFV, species_name, spType, DBH)
  #
  # define tree weights
  tree <- tree %>% mutate(w = case_when(DBH <= 29.9 ~ 10000 / (pi * 7.98^2), DBH > 29.9 ~ 10000 / (pi * 12.62^2)))

  # load list of deciduous and coniferous sp
  deciduousSp <- readRDS('./data/deciduousSp.rds')
  coniferousSp <- readRDS('./data/coniferousSp.rds')
  # define spType
  tree[tree$species_name %in% deciduousSp, 'spType'] <- 'D'
  tree[tree$species_name %in% coniferousSp, 'spType'] <- 'C'
  tree$spType <- as.factor(tree$spType)

  # save tree file
  saveRDS(tree, file = paste0(tempPath, '/treeTemp.rds'))

  ###############################################################
  # load/create data and set extent and resolution
  ###############################################################

  # load mountain mask to remove from case study area
  mountain <- readOGR(dsn = './data/sneznik/GEO', layer = 'Sneznik_area to be deleted', encoding = 'UTF-8', use_iconv = TRUE)
  # enlarge mountain area to get rid of all small patches
  bufferMountain <- buffer(mountain, width = 110, dissolve = T)

  # elevation (m a.s.l.)
  load('./data/sneznik/GEO/topography.Sneznik.rda')
  elevation <- altitude
  names(elevation) <- 'elev'

  # aspect (degrees)
  aspect <- expo_deg
  names(aspect) <- 'aspect'

  # slope (degrees)
  slope <- slope_deg
  names(slope) <- 'slope'

  # get range of Dg and Ba values
  df <- tree %>% group_by(idp) %>%
                 summarise(Dg = sqrt(sum(DBH^2 * w)/sum(w)),
                           BA = sum((pi * (DBH/200)^2) * w))
  #

  # Quadratic diameter (cm)
  dg <- raster('./data/sneznik/GEO/map.DBH.stratified.tif')
  # limit values to range in inventory data
  dg[dg > max(df$Dg)] <- max(df$Dg)
  dg <- resample(dg, elevation)
  names(dg) <- 'Dg'
  # remove mountain area
  dg <- mask(dg, bufferMountain, inverse = TRUE)
  crs(dg) <- crs(mountain)

  # basal area (m2)
  BA <- raster('./data/sneznik/GEO/map.Basal_area.stratified.tif')
  # limit values to range in inventory data
  BA[BA > max(df$BA)] <- max(df$BA)
  BA <- resample(BA, elevation)
  names(BA) <- 'BA'
  # remove mountain area
  BA <- mask(BA, bufferMountain, inverse = TRUE)
  crs(BA) <- crs(mountain)

  # Deciduous proportion
  Dprop <- raster('./data/sneznik/GEO/map.Deciduous_proportion.stratified.tif')
  Dprop <- resample(Dprop, elevation)
  names(Dprop) <- 'Dprop'
  # remove mountain area
  Dprop <- mask(Dprop, bufferMountain, inverse = TRUE)
  crs(Dprop) <- crs(mountain)

  # case study area extent ('park' layer)
  # retrieve from lidar dg map
  parkRaster <- dg
  # exclude area
  parkRaster <- mask(parkRaster, bufferMountain, inverse = TRUE)
  # convert non NA values in 1
  parkRaster[!is.na(parkRaster)] <- as.factor(1)
  names(parkRaster) <- 'park'
  # fill holes in case study area
  # convert raster to SpatRaster
  rasterra <- terra::rast(parkRaster)
  # polygonise
  poly <- terra::as.polygons(rasterra)
  # , trunc=TRUE, dissolve=TRUE, values=TRUE, na.rm=TRUE, extent=FALSE)
  # fill holes
  polyfull <- terra::fillHoles(poly, inverse=FALSE)
  # dissolve smaller polygons into larger ones
  polyagg <- terra::aggregate(polyfull)
  plot(polyagg)
  # convert back to raster
  parkRaster <- terra::rasterize(polyagg, rasterra)
  # convert NA into 0
  isBecomes <- cbind(c(1, NA),
                     c(1, 0))
  parkRaster <- terra::classify(parkRaster, rcl = isBecomes)
  parkRaster <- raster(parkRaster)
  names(parkRaster) <- 'park'

  # create cell ID raster
  ext <- extent(elevation)
  cellID25 <- raster(ext, res = res(elevation))
  cellID25$cellID25 <- c(1:(nrow(cellID25) * ncol(cellID25)))
  cellID25 <- cellID25$cellID25

  ###############################################################
  # save ascii
  ###############################################################

  writeRaster(parkRaster, filename = paste0(landPath, '/parkMask.asc'), format = 'ascii', overwrite = TRUE)
  writeRaster(elevation, filename = paste0(landPath, '/elev.asc'), format = 'ascii', overwrite = TRUE)
  writeRaster(slope, filename = paste0(landPath, '/slope.asc'), format = 'ascii', overwrite = TRUE)
  writeRaster(aspect, filename = paste0(landPath, '/aspect.asc'), format = 'ascii', overwrite = TRUE)
  writeRaster(cellID25, filename = paste0(landPath, '/cellID25.asc'), format = 'ascii', overwrite = TRUE)
  writeRaster(dg, filename = paste0(tempPath, '/dg.grd'), format = 'raster', overwrite = TRUE)
  writeRaster(BA, filename = paste0(tempPath, '/BA.grd'), format = 'raster', overwrite = TRUE)
  writeRaster(Dprop, filename = paste0(tempPath, '/Dprop.grd'), format = 'raster', overwrite = TRUE)

  ###############################################################
  # save data frame
  ###############################################################

  # create raster stack
  rasterStack <- stack(cellID25, parkRaster, elevation, slope,
                       aspect)
  plot(rasterStack)

  # convert into data frame
  envdf <- as.data.frame(rasterStack)

  # save
  envdf$cellID25 <- as.integer(envdf$cellID25)
  envdf$park <- as.integer(envdf$park)
  envdf$elev <- round(envdf$elev, 2)
  envdf$slope <- round(envdf$slope, 2)
  envdf$aspect <- round(envdf$aspect, 2)
  saveRDS(envdf, file = paste0(tempPath, '/envVariablesTemp.rds'))


}
