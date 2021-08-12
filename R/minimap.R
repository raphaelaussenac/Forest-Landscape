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

  # select subset of 1ha cells
  if (landscape == 'bauges'){
    ext <- extent(945000, 945500, 6513400, 6513900) # define random extent
    # make sure to select entire 1ha cells
    cellID100 <- crop(cellID100, ext)
    ext <- extent(cellID100)
  }

  aspect <- crop(aspect, ext)
  cellID25 <- crop(cellID25, ext)
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
  # 3 - transform tree, environment and management data frames
  ###############################################################

  # retrieve cellIDs
  cell100 <- values(cellID100)
  cell25 <- values(cellID25)

  # reduce envVariables df
  df <- read.csv(paste0(landPath, '/envVariables.csv'))
  df <- df[df$cellID100 %in% cell100,]
  # save
  write.csv(df, paste0(miniPath, '/envVariables.csv'), row.names = F)

  # reduce trees75 df
  tree <- read.csv(paste0(landPath, '/trees75.csv'))
  tree <- tree[tree$cellID25 %in% cell25,]
  write.csv(tree, paste0(miniPath, '/trees75.csv'), row.names = F)

  # reduce management df
  manag <- read.csv(paste0(landPath, '/managTable.csv'))
  manag <- manag[manag$cellID100 %in% cell100,]
  write.csv(manag, paste0(miniPath, '/managTable.csv'), row.names = F)

  # surface (ha)
  nrow(df)*625/10000

}
