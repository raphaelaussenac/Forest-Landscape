saveLandscape <- function(landscape){

  ###############################################################
  # initialisation
  ###############################################################

  # load packages
  require(dplyr)
  require(tidyr)
  require(raster)
  require(rgdal)

  # load tree and env data
  tree <- readRDS(paste0(tempPath, '/treeTemp.rds'))
  results <- readRDS(paste0(tempPath, '/trees.rds'))
  envdf <- readRDS(file = paste0(tempPath, '/envVariablesTemp.rds'))

  # load rasters
  cellID25 <- raster(paste0(landPath, '/cellID25.asc'))
  park <- raster(paste0(landPath, '/parkMask.asc'))
  # only consider public forests at milicz
  if(landscape == 'milicz'){
    # load ownership
    own <- readOGR(dsn = './data/milicz/GEO', layer = 'Milicz_forest ownership', encoding = 'UTF-8', use_iconv = TRUE)
  }

  ################################################################################
  # make sure there are exactly the right number of cells in the tree df
  ################################################################################

  # number of cells in the tree data
  nc_tree <- length(unique(results$cellID25))

  # number of cells in the raster data
  nc_raster <- nrow(cellID25) * ncol(cellID25)

  # print error message if cells are duplicated
  if(nc_tree > nc_raster){
    stop("More cells in the tree data than in the raster data. Error could come from the splitting of rasters into strips (in dendro.R), some cells might have been duplicated. Try to split rasters in a different number of strips")
  }
  # print error message if cells are missing
  if(nc_tree < nc_raster){
    stop("Less cells in the tree data than in the raster data. Error could come from the splitting of rasters into strips (in dendro.R), some cells might have been lost. Try to split rasters in a different number of strips")
  }

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
  # create forest column in cell25.csv, forestMask raster and make sure
  # it corresponds to the tree data
  ################################################################################

  # remove cells with no trees in tree data set
  results <- results[!is.na(results$dbh),]

  # remove trees outside park
  results <- left_join(results, envdf %>% dplyr::select(cellID25, park), by = 'cellID25') %>%
             filter(park == 1) %>% dplyr::select(-park)

  # only consider public forests at milicz
  if(landscape == 'milicz'){
    # convert ownership into raster
    own <- rasterize(own[,1], cellID25, background = 10)
    names(own) <- 'public'
    # convert values
    isBecomes <- cbind(c(1, 2, 10),
                       c(1, 0, 0))
    own <- reclassify(own, rcl = isBecomes)
    names(own) <- 'public'
    # stack with cellID25
    own <- stack(cellID25, own)
    # convert into dataframe
    own <- as.data.frame(own)
    # remove trees outside public forests
    results <- left_join(results, own, by = 'cellID25') %>% filter(public == 1) %>% dplyr::select(-public)
  }

  # add cellID100 to envdf
  cellID <- stack(cellID25, cellID100)
  cellID <- as.data.frame(values(cellID))
  envdf <- full_join(envdf, cellID, by = 'cellID25')

  # add forest cells to envdf
  listForestCells <- unique(results$cellID25)
  envdf$forest <- NA
  envdf[envdf$cellID25 %in% listForestCells, 'forest'] <- 1
  envdf[!(envdf$cellID25 %in% listForestCells), 'forest'] <- 0

  # reduce table size in memory
  envdf$cellID100 <- as.integer(envdf$cellID100)
  envdf$forest <- as.integer(envdf$forest)
  # sort colnames
  colOrd <- c('cellID25','cellID100','park', 'forest', 'elev','slope','aspect')
  if(landscape == 'bauges'){
    colOrd <- c(colOrd,'swhc','pH','GRECO','SIQpet','SIFsyl','SIAalb','SIPabi')
  }
  # format
  envdf <- envdf %>% dplyr::select(all_of(colOrd)) %>% arrange(cellID25, cellID100) %>%
                     mutate(across(c('elev','slope','aspect'), round, 4))
  # save
  write.csv(envdf, file = paste0(landPath, '/cell25.csv'), row.names = FALSE)

  # save forestMask raster
  cellID25$forest <- envdf$forest
  forestMask <- cellID25$forest
  writeRaster(forestMask, filename = paste0(landPath, '/forestMask.asc'), format = 'ascii', overwrite = TRUE)

  # reduce table size in memory
  results$cellID25 <- as.integer(results$cellID25)
  results$n <- as.integer(results$n)
  results$dbh <- round(results$dbh, 2)

  # save
  saveRDS(results[, c('cellID25', 'sp', 'n', 'dbh', 'wlid')], file = paste0(tempPath, '/trees75.rds'))

}
