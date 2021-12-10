prepMilicz <- function(){

  ###############################################################
  # initialisation
  ###############################################################

  # load packages
  require(rgdal)
  require(raster)
  # rasterize is faster with velox but package is
  # not yet compatible with R 4.0.4
  # require(velox)
  # require(gdalUtils)
  require(sf)
  require(stringr)
  require(dplyr)

  ###############################################################
  # load tree data
  ###############################################################

  # load tree inventory data
  tree <- read.csv('./data/milicz/inventory/Milicz_trees.csv', sep = ';')

  # function to replace commas by dots and transform character vectors
  # into numeric vectors
  comTodot <- function(charVect){
    charVect <- as.numeric(str_replace_all(charVect, ',', '\\.'))
    return(charVect)
  }

  tree$DBH <- comTodot(tree$DBH)
  tree$height <- comTodot(tree$height)
  tree$species <- as.factor(tree$species)
  tree$species <- recode_factor(tree$species, 'Malus silvestris' = 'Malus sylvestris',
                                                  'Quercus undefined' = 'Quercus robur',
                                                  'Rhamnus frangula' = 'Frangula alnus',
                                                  'Ulmus' = 'Ulmus minor')
  #

  # add missing columns remove useless columns
  tree <- tree %>% dplyr::select(plot_id, DBH, species) %>%
                   rename(idp = plot_id, species_name = species) %>%
                   mutate(idp = as.factor(idp), w = NA, spType = NA, CODE_TFV = as.factor(99)) %>%
                   relocate(idp, w, CODE_TFV, species_name, spType, DBH)
  #
  # define tree weights
  # plot radius = 12.62 m --> pi*12.62^2 = 500.3439 m2
  tree$w <- 10000 / 500.3439

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

  # elevation (m a.s.l.)
  load('./data/milicz/GEO/topography.Milicz.rda')
  elevation <- altitude
  names(elevation) <- 'elev'

  # aspect (degrees)
  aspect <- expo_deg
  names(aspect) <- 'aspect'

  # slope (degrees)
  slope <- slope_deg
  names(slope) <- 'slope'

  # case study area extent ('park' layer)
  # retrieve from lidar slope map
  parkRaster <- slope
  # convert non NA values in 1
  parkRaster[!is.na(parkRaster)] <- 1
  # convert NA into 0
  isBecomes <- cbind(c(1, NA),
                     c(1, 0))
  parkRaster <- reclassify(parkRaster, rcl = isBecomes)
  names(parkRaster) <- 'park'

  # get range of Dg and Ba values
  df <- tree %>% group_by(idp) %>%
                 summarise(Dg = sqrt(sum(DBH^2 * w)/sum(w)),
                           BA = sum((pi * (DBH/200)^2) * w))
  #

  # Quadratic diameter (cm)
  dg <- raster('./data/milicz/GEO/map.DBH.stratified.tif')
  # limit values to range in inventory data
  dg[dg > max(df$Dg)] <- max(df$Dg)
  crs(dg) <- crs(parkRaster)
  dg <- resample(dg, elevation)
  names(dg) <- 'Dg'

  # basal area (m2)
  BA <- raster('./data/milicz/GEO/map.BA.stratified.tif')
  # limit values to range in inventory data
  BA[BA > max(df$BA)] <- max(df$BA)
  crs(BA) <- crs(parkRaster)
  BA <- resample(BA, elevation)
  names(BA) <- 'BA'

  # Deciduous proportion
  Dprop <- raster('./data/milicz/GEO/map.DP.stratified.tif')
  crs(Dprop) <- crs(parkRaster)
  Dprop <- resample(Dprop, elevation)
  names(Dprop) <- 'Dprop'

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
