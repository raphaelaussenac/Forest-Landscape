saveLandscape <- function(landscape){

  ###############################################################
  # initialisation
  ###############################################################

  # load packages
  require(dplyr)
  require(tidyr)
  require(raster)

  # load tree and env data
  tree <- readRDS(paste0(tempPath, '/treeTemp.rds'))
  results <- readRDS(paste0(tempPath, '/trees.rds'))
  envdf <- readRDS(file = paste0(tempPath, '/envVariablesTemp.rds'))

  # load rasters
  cellID25 <- raster(paste0(landPath, '/cellID25.asc'))
  compo <- raster(paste0(tempPath, '/compoID.asc'))
  park <- raster(paste0(landPath, '/parkMask.asc'))

  ################################################################################
  # remove trees with dbh < min dbh of tree-level inventory data
  ################################################################################

  # define min dbh of tree-level inventory data
  minDBH <- min(tree$DBH)

  # remove trees < minDBH
  # if all trees are < minDBH then replace by a line of NA

  # first count total number of lines per cell, number of line with dbh < minDBH
  # and number of lines with dbh >= minDBH
  treeCnt <- results %>% filter(!is.na(sp)) %>% mutate(n = 1, size = ifelse(dbh >= minDBH, 'big', 'small')) %>%
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
  results <- results %>% filter(!( cellID25 %in% removelist$cellID25 & dbh <minDBH)) %>% dplyr::select(-i)

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
  writeRaster(cellID100, filename = paste0(landPath, '/cellID100.asc'), format = 'ascii', overwrite = TRUE)
  # change resolution back to cellID25 resolution in order to crop and then
  # stack together
  cellID100 <- disaggregate(cellID100, fact = resFact)
  # reduce extent of cellID100 to minimum common extent with cellID25
  cellID100 <- crop(cellID100, extent(cellID25))


  ################################################################################
  # create forest extent raster:
  ################################################################################

  # wherever there's a compoID inside the case study area (= park)
  compo[!is.na(compo)] <- 1
  # convert NA into 0
  isBecomes <- cbind(c(1, NA),
                     c(1, 0))
  compo <- reclassify(compo, rcl = isBecomes)
  # union of park and compo
  compoUpark <- park + compo
  # convert values
  isBecomes <- cbind(c(0, 1, 2),
                     c(0, 0, 1))
  forest <- reclassify(compoUpark, rcl = isBecomes)
  names(forest) <- 'forest'
  # add forest to cell data
  cellID <- stack(cellID25, cellID100, forest)
  plot(cellID)
  writeRaster(cellID$forest, filename = paste0(landPath, '/forestMask.asc'), format = 'ascii', overwrite = TRUE)

  # add cellID100 and forest extent to env data
  cellIDforest <- as.data.frame(values(cellID))
  envdf <- merge(envdf, cellIDforest, by = 'cellID25')

  ################################################################################
  # save
  ################################################################################

  # reduce table size in memory
  envdf$cellID100 <- as.integer(envdf$cellID100)
  envdf$forest <- as.integer(envdf$forest)
  # sort colnames
  colOrd <- c('cellID25','cellID100','park', 'forest', 'elev','slope','aspect')
  if(landscape == 'bauges'){
    colOrd <- c(colOrd,'swhc','pH','GRECO','SIQpet','SIFsyl','SIAalb','SIPabi')
  }
  # save
  envdf <- envdf %>% dplyr::select(all_of(colOrd)) %>% arrange(cellID25, cellID100) %>%
                     mutate(across(c('elev','slope','aspect'), round, 4))
  write.csv(envdf, file = paste0(landPath, '/cell25.csv'), row.names = FALSE)

  # remove cells with no trees
  results <- results[!is.na(results$dbh),]
  # remove trees not in forest cells
  forestCells <- envdf[envdf$forest == 1, 'cellID25']
  results <- results %>% filter(cellID25 %in% forestCells)
  # reduce table size in memory
  results$cellID25 <- as.integer(results$cellID25)
  results$n <- as.integer(results$n)
  results$dbh <- round(results$dbh, 2)

  # save
  saveRDS(results[, c('cellID25', 'sp', 'n', 'dbh', 'wlid')], file = paste0(tempPath, './trees75.rds'))

}
