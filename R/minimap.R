minimap <- function(landscape){

  ###############################################################
  # initialisation
  ###############################################################

  # load packages
  require(raster)
  require(dplyr)

  ###############################################################
  # 1 - open all rasters
  ###############################################################

  aspect <- raster(paste0(landPath, '/aspect.asc'))
  cellID25 <- raster(paste0(landPath, '/cellID25.asc'))
  cellID100 <- raster(paste0(landPath, '/cellID100.asc'))
  elev <- raster(paste0(landPath, '/elev.asc'))
  forest <- raster(paste0(landPath, '/forestMask.asc'))
  park <- raster(paste0(landPath, '/parkMask.asc'))
  pH <- raster(paste0(landPath, '/pH.asc'))
  slope <- raster(paste0(landPath, '/slope.asc'))
  swhc <- raster(paste0(landPath, '/swhc.asc'))

  ###############################################################
  # 2 - crop on smaller extent
  ###############################################################

  if (landscape == 'bauges'){
    ext <- extent(945000, 945500, 6513400, 6513900)
  }

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

  writeRaster(park, filename = paste0(miniPath, '/parkMask.asc'), format = 'ascii', overwrite = TRUE)
  writeRaster(forest, filename = paste0(miniPath, '/forestMask.asc'), format = 'ascii', overwrite = TRUE)
  writeRaster(elev, filename = paste0(miniPath, '/elev.asc'), format = 'ascii', overwrite = TRUE)
  writeRaster(slope, filename = paste0(miniPath, '/slope.asc'), format = 'ascii', overwrite = TRUE)
  writeRaster(aspect, filename = paste0(miniPath, '/aspect.asc'), format = 'ascii', overwrite = TRUE)
  writeRaster(swhc, filename = paste0(miniPath, '/swhc.asc'), format = 'ascii', overwrite = TRUE)
  writeRaster(pH, filename = paste0(miniPath, '/pH.asc'), format = 'ascii', overwrite = TRUE)
  writeRaster(cellID25, filename = paste0(miniPath, '/cellID25.asc'), format = 'ascii', overwrite = TRUE)
  writeRaster(cellID100, filename = paste0(miniPath, '/cellID100.asc'), format = 'ascii', overwrite = TRUE)

  ###############################################################
  # 3 - transform tree and environmental data frames
  ###############################################################

  # retrieve cellID
  cellID <- values(cellID25)

  # reduce envVariables df
  df <- read.csv(paste0(landPath, '/envVariables.csv'))
  df <- df[df$cellID25 %in% cellID,]
  # calculate nb of forest cells per ha for the mini landscape
  # may be different from nb of forest cells in envVariables.csv
  df <- df %>% group_by(cellID100) %>% mutate(forestCellsPerHa = sum(forest))
  # plot hist of nb of forest cells per ha
  test <- df %>% group_by(cellID100) %>% summarise(forestCellsPerHa = unique(forestCellsPerHa))
  hist(test$forestCellsPerHa, breaks = c(0:16))
  # save
  write.csv(df, paste0(miniPath, '/envVariables.csv'), row.names = F)

  # surface (ha)
  nrow(df)*625/10000

  # reduce trees75 df
  df <- read.csv(paste0(landPath, '/trees75.csv'))
  df <- df[df$cellID25 %in% cellID,]
  write.csv(df, paste0(miniPath, '/trees75.csv'), row.names = F)


}
