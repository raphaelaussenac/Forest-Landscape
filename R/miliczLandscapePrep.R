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
  require(gdalUtils)
  require(sf)

  ###############################################################
  # load/create data and set extent and resolution
  ###############################################################

  # forest
  forest <- readOGR(dsn = './data/milicz/GEO', layer = 'Milicz_forest_stands', encoding = 'UTF-8', use_iconv = TRUE)

  # elevation (m a.s.l.)
  load('./data/milicz/GEO/altitude.rda')
  elevation <- altitude
  names(elevation) <- 'elev'
  # set projection
  crs(elevation) <- crs(forest)

  # slope (degree)
  slope <- terrain(elevation, opt = 'slope', unit = 'degrees', neighbors = 8)

  # aspect (degree)
  aspect <- terrain(elevation, opt = 'aspect', unit = 'degrees', neighbors = 8)

  # convert forest into raster and set extent + resolution
  ext <- extent(elevation)
  r <- raster(ext, res = res(elevation))
  forest$ID <- as.factor(forest$stand_id)
  forestRaster <- rasterize(forest, r, field = 'ID')
  # set projection
  crs(forestRaster) <- crs(forest)
  # convert NA into 0
  isBecomes <- cbind(c(1:nrow(forest), NA),
                     c(rep(1, nrow(forest)), 0))
  forestRaster <- reclassify(forestRaster, rcl = isBecomes)
  names(forestRaster) <- 'forest'

  # Quadratic diameter (cm) [0, 80]
  dg <- raster('./data/milicz/GEO/map.DBH.stratified.tif')
  crs(dg) <- crs(forest)
  dg <- resample(dg, elevation)

  # basal area (m2) [0, 120]
  BA <- raster('./data/milicz/GEO/map.BA.stratified.tif')
  crs(BA) <- crs(forest)
  BA <- resample(BA, elevation)

  # Deciduous proportion (% of total BA)
  Dprop <- raster('./data/milicz/GEO/map.DP.stratified.tif')
  crs(Dprop) <- crs(forest)
  Dprop <- resample(Dprop, elevation)

  # create cell ID raster
  ext <- extent(elevation)
  cellID25 <- raster(ext, res = res(elevation))
  cellID25$cellID25 <- c(1:(nrow(cellID25) * ncol(cellID25)))
  cellID25 <- cellID25$cellID25

  ###############################################################
  # save ascii
  ###############################################################

  writeRaster(forestRaster, filename = paste0(landPath, '/forestMask.asc'), format = 'ascii', overwrite = TRUE)
  writeRaster(elevation, filename = paste0(landPath, '/elev.asc'), format = 'ascii', overwrite = TRUE)
  writeRaster(slope, filename = paste0(landPath, '/slope.asc'), format = 'ascii', overwrite = TRUE)
  writeRaster(aspect, filename = paste0(landPath, '/aspect.asc'), format = 'ascii', overwrite = TRUE)
  writeRaster(cellID25, filename = paste0(landPath, '/cellID25.asc'), format = 'ascii', overwrite = TRUE)
  writeRaster(dg, filename = paste0(tempPath, '/dg.asc'), format = 'ascii', overwrite = TRUE)
  writeRaster(BA, filename = paste0(tempPath, '/BA.asc'), format = 'ascii', overwrite = TRUE)
  writeRaster(Dprop, filename = paste0(tempPath, '/Dprop.asc'), format = 'ascii', overwrite = TRUE)

  ###############################################################
  # save data frame
  ###############################################################

  # create raster stack
  rasterStack <- stack(cellID25, forestRaster, elevation, slope,
                       aspect)
  plot(rasterStack)

  # convert into data frame
  envdf <- as.data.frame(rasterStack)

  # save
  envdf$cellID25 <- as.integer(envdf$cellID25)
  envdf$forest <- as.integer(envdf$forest)
  envdf$elev <- round(envdf$elev, 2)
  envdf$slope <- round(envdf$slope, 2)
  envdf$aspect <- round(envdf$aspect, 2)
  saveRDS(envdf, file = paste0(tempPath, '/envVariablesTemp.rds'))

}
