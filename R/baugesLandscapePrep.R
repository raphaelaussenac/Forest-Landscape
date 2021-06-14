prepBauges <- function(){

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

  # park
  park <- readOGR(dsn = './data/bauges/GEO', layer = 'park', encoding = 'UTF-8', use_iconv = TRUE)

  # forest
  forest <- readOGR(dsn = './data/bauges/GEO', layer = 'BD_Foret_V2_PNRfilled_Foret_2014', encoding = 'UTF-8', use_iconv = TRUE)

  # elevation (m a.s.l.)
  elevation <- raster('./data/bauges/GEO/MNT_all_5m.tif')
  # set projection
  crs(elevation) <- crs(park)
  # set resolution to 25*25 instead of 5*5
  elevation <- aggregate(elevation, fact = 5)
  names(elevation) <- 'elev'

  # slope (degree)
  slope <- terrain(elevation, opt = 'slope', unit = 'degrees', neighbors = 8)

  # aspect (degree)
  aspect <- terrain(elevation, opt = 'aspect', unit = 'degrees', neighbors = 8)

  # convert park into raster and set extent + resolution
  ext <- extent(elevation)
  r <- raster(ext, res=res(elevation))
  park$ID <- as.factor(park$ID)
  parkRaster <- rasterize(park, r, field = 'ID')
  # set projection
  crs(parkRaster) <- crs(park)
  # convert NA into 0
  isBecomes <- cbind(c(1, NA),
                     c(1, 0))
  parkRaster <- reclassify(parkRaster, rcl = isBecomes)
  names(parkRaster) <- 'park'

  # convert forest into raster and set extent + resolution
  forest$ID <- as.factor(forest$ID)
  forestRaster <- rasterize(forest, r, field = 'ID')
  # set projection
  crs(forestRaster) <- crs(forest)
  # convert NA into 0
  isBecomes <- cbind(c(1:nrow(forest), NA),
                     c(rep(1, nrow(forest)), 0))
  forestRaster <- reclassify(forestRaster, rcl = isBecomes)
  names(forestRaster) <- 'forest'

  # swhc (cm)
  swhc <- raster('./data/bauges/GEO/rum_500_v2009.tif')
  swhc <- resample(swhc, elevation)
  swhc <- swhc/10 # convert into cm
  names(swhc) <- 'swhc'

  # pH
  pH <- raster('./data/bauges/GEO/ph_2008.tif')
  pH <- resample(pH, elevation)
  names(pH) <- 'pH'

  # Quadratic diameter (cm) [0, 80]
  dg <- raster('./data/bauges/GEO/rastDg75error.clean.tif')
  crs(dg) <- crs(park)
  dg <- resample(dg, elevation)

  # basal area (m2) [0, 120]
  BA <- raster('./data/bauges/GEO/rastG75error.clean.tif')
  crs(BA) <- crs(park)
  BA <- resample(BA, elevation)

  # # Quadratic diameter (cm) with all extreme values
  # dgExtrem <- raster('./data/bauges/GEO/rastDg75error.forest.tif')
  # crs(dgExtrem) <- crs(park)
  # dgExtrem <- resample(dgExtrem, elevation)
  #
  # # basal area (m2) with all extreme values
  # BAExtrem <- raster('./data/bauges/GEO/rastG75error.forest.tif')
  # crs(BAExtrem) <- crs(park)
  # BAExtrem <- resample(BAExtrem, elevation)

  # N
  N <- raster('./data/bauges/GEO/N_pred.tif')
  crs(N) <- crs(park)
  N <- resample(N, elevation)

  # large tree (BA or proportion of total BA?)
  # LTBA <- raster('./data/bauges/GEO/GGB_pred.tif')
  # crs(LTBA) <- crs(park)
  # LTBA <- resample(LTBA, elevation)

  # Deciduous proportion (% of total BA)
  Dprop <- raster('./data/bauges/GEO/propGR_ONF_25filled.tif')
  Dprop <- 100 - Dprop
  crs(Dprop) <- crs(park)
  Dprop <- resample(Dprop, elevation)

  # create cell ID raster
  ext <- extent(elevation)
  cellID25 <- raster(ext, res=res(elevation))
  cellID25$cellID25 <- c(1:(nrow(cellID25) * ncol(cellID25)))
  cellID25 <- cellID25$cellID25

  # greco
  greco <- readOGR(dsn = './data/bauges/GEO', layer = 'greco_l93', encoding = 'UTF-8', use_iconv = TRUE)
  # convert into a raster
  greco$CODEGRECO <- as.factor(greco$CODEGRECO)
  grecoRaster <- rasterize(greco, cellID25, field='CODEGRECO')
  crs(grecoRaster) <- crs(park)
  names(grecoRaster) <- 'GRECO'

  ###############################################################
  # save ascii
  ###############################################################

  writeRaster(parkRaster, filename = paste0(landPath, '/parkMask.asc'), format = 'ascii', overwrite = TRUE)
  writeRaster(forestRaster, filename = paste0(landPath, '/forestMask.asc'), format = 'ascii', overwrite = TRUE)
  writeRaster(elevation, filename = paste0(landPath, '/elev.asc'), format = 'ascii', overwrite = TRUE)
  writeRaster(slope, filename = paste0(landPath, '/slope.asc'), format = 'ascii', overwrite = TRUE)
  writeRaster(aspect, filename = paste0(landPath, '/aspect.asc'), format = 'ascii', overwrite = TRUE)
  writeRaster(swhc, filename = paste0(landPath, '/swhc.asc'), format = 'ascii', overwrite = TRUE)
  writeRaster(pH, filename = paste0(landPath, '/pH.asc'), format = 'ascii', overwrite = TRUE)
  writeRaster(cellID25, filename = paste0(landPath, '/cellID25.asc'), format = 'ascii', overwrite = TRUE)
  writeRaster(dg, filename = paste0(tempPath, '/dg.asc'), format = 'ascii', overwrite = TRUE)
  writeRaster(BA, filename = paste0(tempPath, '/BA.asc'), format = 'ascii', overwrite = TRUE)
  writeRaster(N, filename = paste0(tempPath, '/N.asc'), format = 'ascii', overwrite = TRUE)
  writeRaster(Dprop, filename = paste0(tempPath, '/Dprop.asc'), format = 'ascii', overwrite = TRUE)
  writeRaster(grecoRaster, filename = paste0(tempPath, '/greco.asc'), format = 'ascii', overwrite = TRUE)

  ###############################################################
  # save data frame
  ###############################################################

  # create raster stack
  rasterStack <- stack(cellID25, parkRaster, forestRaster, elevation, slope,
                       aspect, swhc, pH, grecoRaster)
  plot(rasterStack)

  # convert into data frame
  envdf <- as.data.frame(rasterStack)

  # rename back GRECO regions
  envdf[envdf$GRECO == 4, 'GRECO'] <- 'C'
  envdf[envdf$GRECO == 6, 'GRECO'] <- 'E'
  envdf[envdf$GRECO == 9, 'GRECO'] <- 'H'

  # predict salem SI (site index)
  envdf <- salemSI(envdf)

  # save
  envdf$cellID25 <- as.integer(envdf$cellID25)
  envdf$park <- as.integer(envdf$park)
  envdf$forest <- as.integer(envdf$forest)
  envdf$elev <- round(envdf$elev, 2)
  envdf$slope <- round(envdf$slope, 2)
  envdf$aspect <- round(envdf$aspect, 2)
  envdf$swhc <- round(envdf$swhc, 2)
  envdf$pH <- round(envdf$pH, 2)
  envdf$GRECO <- as.factor(envdf$GRECO)
  envdf$SIQpet <- round(envdf$SIQpet, 2)
  envdf$SIFsyl <- round(envdf$SIFsyl, 2)
  envdf$SIAalb <- round(envdf$SIAalb, 2)
  envdf$SIPabi <- round(envdf$SIPabi, 2)
  saveRDS(envdf, file = paste0(tempPath, '/envVariablesTemp.rds'))


  ###############################################################
  # save SI maps and distributions
  ###############################################################

  # pred SI map
  rasterStack$SIQpet <- envdf$SIQpet
  rasterStack$SIFsyl <- envdf$SIFsyl
  rasterStack$SIAalb <- envdf$SIAalb
  rasterStack$SIPabi <- envdf$SIPabi

  pdf(paste0(evalPath, '/salemSIpredMap.pdf'), width = 10, height = 10)
  plot(rasterStack[[10:13]])
  dev.off()

  # pred SI distribution
  pdf(paste0(evalPath, '/salemSIpredDist.pdf'), width = 10, height = 10)
  par(mfrow = c(2,2))
  hist(envdf$SIQpet, main = 'Q. petraea')
  hist(envdf$SIFsyl, main = 'F. sylvatica')
  hist(envdf$SIAalb, main = 'A. alba')
  hist(envdf$SIPabi, main = 'P. abies')
  dev.off()

}
