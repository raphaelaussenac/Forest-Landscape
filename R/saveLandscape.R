saveLandscape <- function(landscape){

  ###############################################################
  # initialisation
  ###############################################################

  # load packages
  require(dplyr)
  require(tidyr)
  require(raster)

  # load tree data
  results <- readRDS(paste0(tempPath, '/trees.rds'))

  # load cellID25 raster
  cellID25 <- raster(paste0(landPath, '/cellID25.asc'))

  ################################################################################
  # remove trees <7.5 cm
  ################################################################################

  # remove trees <7.5cm in cells where there are some trees left
  # otherwise, if all trees are <7.5 then replace by a line of NA

  # first count total number of lines per cell, number of line with dbh <7.5
  # and number of lines with dbh >= 7.5
  treeCnt <- results %>% filter(!is.na(sp)) %>% mutate(n = 1, size = ifelse(dbh >= 7.5, 'big', 'small')) %>%
                     group_by(cellID25) %>% count(size, wt = n) %>% pivot_wider(id_cols = cellID25, names_from = size, values_from = n) %>%
                     mutate(nbtot = sum(big, small, na.rm = TRUE))
  # list of cells where small trees can (simply) be removed
  removelist <- treeCnt %>% filter(small > 0 & !is.na(small) & big > 0 & !is.na(big))
  removelist <- as.list(removelist[, 'cellID25'])

  # list of cells where small trees make up the whole stand
  # and cannot be removed --> create a NA line instead
  setToNAlist <- treeCnt %>% filter(small > 0 & !is.na(small) & is.na(big))
  setToNAlist <- as.list(setToNAlist[, 'cellID25'])

  # remove trees
  results <- results %>% filter(!( cellID25 %in% removelist$cellID25 & dbh <7.5)) %>% dplyr::select(-i)

  # set to NA
  # for that we remove the cells and add new NA lines
  results <- results %>% filter(! (cellID25 %in% setToNAlist$cellID25))
  # NAdf <- data.frame(sp = NA, wlid = NA, n = NA, dbh = NA, i = NA, cellID25 = setToNAlist$cellID25)
  NAdf <- data.frame(cellID25 = setToNAlist$cellID25, sp = NA, n = NA, dbh = NA, wlid = NA)
  results <- rbind(results, NAdf)
  results <- results %>% arrange(cellID25)


  ################################################################################
  # aggregate data at 100*100 resolution
  ################################################################################

  # create second cellID raster with different resolution (25*25 -> 100*100)
  resFact <- 4 # multiply resolution by 4
  cellID100 <- aggregate(cellID25, fact = resFact)
  names(cellID100) <- 'cellID100'
  # assign unique cellID
  values(cellID100) <- c(1:(nrow(cellID100) * ncol(cellID100)))
  # save raster cellID100 at 100*100 resolution
  writeRaster(cellID100, filename = paste0(landPath, './cellID100.asc'), format = 'ascii', overwrite = TRUE)
  # change resolution back to cellID25 resolution in order to crop and then
  # stack together
  cellID100 <- disaggregate(cellID100, fact = resFact)
  # reduce extent of cellID100 to minimum common extent with cellID25
  cellID100 <- crop(cellID100, extent(cellID25))
  # stack rasters
  cellID <- stack(cellID25, cellID100)
  # cellID correspondence dataframe
  cellIDdf <- as.data.frame(values(cellID))


  ################################################################################
  # save
  ################################################################################

  # create forest extent raster
  # create df of unique cellID with total number of trees N
  df <- results %>% group_by(cellID25) %>% summarise(N = sum(n))
  # define the forest variable based on N
  df <- df %>% mutate(forest = if_else(is.na(N), 0, 1))
  # add forest to cellID raster
  cellID$forest <- df$forest
  plot(cellID$forest)
  # save raster
  writeRaster(cellID$forest, filename = paste0(landPath, './forestMask.asc'), format = 'ascii', overwrite = TRUE)

  # add cellID100 to environmental data
  envdf <- readRDS(file = paste0(tempPath, '/envVariablesTemp.rds'))
  cellIDforest <- as.data.frame(values(cellID))
  envdf <- merge(envdf, cellIDforest, by = 'cellID25', all.y = F)
  if(landscape != 'bauges'){
    envdf <- envdf %>% dplyr::select(-forest.y) %>% rename(forest = forest.x) %>%
                       mutate(park = NA)
  }
  # reduce table size in memory
  envdf$cellID100 <- as.integer(envdf$cellID100)
  # sort colnames
  colOrd <- c('cellID25','cellID100','park', 'forest', 'elev','slope','aspect')
  if(landscape == 'bauges'){
    colOrd <- c(colOrd,'swhc','pH','GRECO','SIQpet','SIFsyl','SIAalb','SIPabi')
  }
  # save
  write.csv(envdf[, colOrd], file = paste0(landPath, './cell25.csv'), row.names = FALSE)

  # remove cells with no trees
  results <- results[!is.na(results$dbh),]
  # reduce table size in memory
  results$cellID25 <- as.integer(results$cellID25)
  results$n <- as.integer(results$n)
  results$dbh <- round(results$dbh, 2)

  # save
  saveRDS(results[, c('cellID25', 'sp', 'n', 'dbh', 'wlid')], file = paste0(tempPath, './trees75.rds'))

}
