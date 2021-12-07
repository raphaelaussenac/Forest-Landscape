dendro <- function(landscape){

  ###############################################################
  # initialisation
  ###############################################################

  # load packages
  require(dplyr)
  require(raster)
  require(doParallel)
  require(tidyr)

  # load composition ID
  compo <- raster(paste0(tempPath, '/compoID.asc'))

  # load lidar data
  Dg <- raster(paste0(tempPath, '/dg.asc'))
  Dg <- setMinMax(Dg)
  BA <- raster(paste0(tempPath, '/BA.asc'))
  BA <- setMinMax(BA)
  Dprop <- raster(paste0(tempPath, '/Dprop.asc'))
  Dprop <- setMinMax(Dprop)
  cellID25 <- raster(paste0(landPath, '/cellID25.asc'))

  # create raster stack
  rasterStack <- stack(compo, Dg, BA, Dprop, cellID25)

  # load NFI tree data
  tree <- readRDS(paste0(tempPath, '/treeTemp.rds'))

  ###############################################################
  # calculate N, Dg, BA, Dprop in NFI plots
  ###############################################################

  # Calculate N, Dg and BA for the whole plote, for deciduous and coniferous and
  # for ech species
  tree <- tree %>% mutate(BAtree =  (pi * (DBH/200)^2) * w) %>%
                   group_by(idp) %>% mutate(BAtot = sum((pi * (DBH/200)^2) * w),
                                            Dgtot = sqrt(sum(DBH^2 * w)/sum(w)),
                                            Ntot = sum(w)) %>%
                   group_by(idp, spType) %>% mutate(BAdc = sum((pi * (DBH/200)^2) * w),
                                                    Dgdc = sqrt(sum(DBH^2 * w)/sum(w)),
                                                    Ndc = sum(w)) %>%
                   group_by(idp, spType, species_name) %>% mutate(BAsp = sum((pi * (DBH/200)^2) * w),
                                                                  Dgsp = sqrt(sum(DBH^2 * w)/sum(w))) %>%
                   group_by(idp, spType) %>% mutate(spPropdc = BAsp/BAdc) %>% ungroup() %>%
                   mutate(treePropsp = BAtree / BAsp)
  #
  # summarise at sp level
  NFIsp <- tree %>% group_by(idp, species_name, BAsp, Dgsp, spType, spPropdc) %>% summarise()

  ###############################################################
  #
  ###############################################################

  # Split data into 4 separate stack rasters
  xmin <- extent(rasterStack)[1]
  xmax <- extent(rasterStack)[2]
  ystep <- round((extent(rasterStack)[4] - extent(rasterStack)[3]) / 4)
  ymin <- c(extent(rasterStack)[3], extent(rasterStack)[3] + ystep, extent(rasterStack)[3] + ystep * 2 , extent(rasterStack)[3] + ystep *3)
  ymax <- c(extent(rasterStack)[3] + ystep, extent(rasterStack)[3] + ystep * 2 , extent(rasterStack)[3] + ystep * 3, extent(rasterStack)[4])

  rast1 <- crop(rasterStack, c(xmin, xmax, ymin[4], ymax[4]))
  rast2 <- crop(rasterStack, c(xmin, xmax, ymin[3], ymax[3]))
  rast3 <- crop(rasterStack, c(xmin, xmax, ymin[2], ymax[2]))
  rast4 <- crop(rasterStack, c(xmin, xmax, ymin[1], ymax[1]))

  # i <- 37150 # milicz
  # cell <- rast3[i]
  # cell

  # function to assign n trees and their dbh to each cell
  assignDendro <- function(cell, i, tree, NFIsp){

    # retrieve id composition dg and ba (LIDAR)
    id <- cell[1]
    dgtot <- cell[2]
    batot <- cell[3]
    dprop <- cell[4]
    cellID25 <- cell[5]

    if(sum(is.na(cell)) == 0 & dgtot > 0 & batot > 0){
      # deciduous - coniferous level -------------------------------------------------
      # calculate deciduous and coniferous ba (LIDAR)
      bad <- batot * dprop
      bac <- batot * (1 - dprop)

      # if there is only deciduous or coniferous species in the NFI plot
      # set bad and bac to 0 and 100 accordingly
      if(nrow(tree[tree$idp == id & tree$spType == 'D',]) == 0){
        bad <- batot * 0
        bac <- batot * (1 - 0)
      }
      if(nrow(tree[tree$idp == id & tree$spType == 'C',]) == 0){
        bad <- batot * 1
        bac <- batot * (1 - 1)
      }

      # retrieve dg deciduous and coniferous (NFI)
      dgd <- as.numeric(unique(tree[tree$idp == id & tree$spType == 'D', 'Dgdc']))
      dgc <- as.numeric(unique(tree[tree$idp == id & tree$spType == 'C', 'Dgdc']))

      # calculate alpha correction coef for deciduous and coniferous
      alpha <- dgtot * sqrt( sum(bad/dgd^2, bac/dgc^2, na.rm = TRUE) / batot)

      # sp level ---------------------------------------------------------------------

      # calculate species ba (LIDAR) from NFI sp proportion
      NFIplot <- NFIsp[NFIsp$idp == id,]
      NFIplot[NFIplot$spType == 'D', 'BAdclid'] <- bad
      NFIplot[NFIplot$spType == 'C', 'BAdclid'] <- bac
      NFIplot$BAsplid <- NFIplot$BAdclid * NFIplot$spPropdc

      # tree level -------------------------------------------------------------------

      # calculate tree ba (LIDAR) from NFI tree ba proportion
      dbh <- tree[tree$idp == id,]
      dbh <- merge(dbh, NFIplot[, c('species_name', 'BAsplid')], by = 'species_name')
      dbh$BAtreelid <- dbh$BAsplid * dbh$treePropsp

      # assign a DBH (LIDAR) to all trees
      dbh$DBHlid <- dbh$DBH * alpha

      # assign weight to each tree
      dbh$wlid <- 40000 / pi * dbh$BAtreelid / dbh$DBHlid^2

      # add i and cellID now ------
      dbh$i <- i
      dbh$cellID25 <- cellID25
      # ---------------------------

      # calculate round(weight) for a 25*25m pixel --------------------------------------------- method 1
      # and set min weight to 1 or 0
      # dbh$w25m <- apply(as.data.frame(dbh[, 'wlid']), 1, function(x) max(1, round(x/16)))

      # # calculate weight for a 25*25m pixel (with stochastic process for w <= 0.5) ----------- method 2
      # dbh$w25m <- dbh$wlid / 16
      # # if weight > 0.5 --> round(weight)
      # dbh[dbh$w25m > 0.5, 'w25m'] <- round(dbh[dbh$w25m > 0.5, 'w25m'])
      # # if weight <= 0.5 --> set weight to 0 or 1 depending on weight
      # # (stochastic process)
      # dbh$random <- runif(nrow(dbh), min = 0, max = 0.5)
      # dbh$randSmallerThanw25m <- dbh$random < dbh$w25m
      # dbh[dbh$w25m <= 0.5 & dbh$randSmallerThanw25m == TRUE, 'w25m'] <- 1
      # dbh[dbh$w25m <= 0.5 & dbh$randSmallerThanw25m == FALSE, 'w25m'] <- 0

      # # calculate weight for a 25*25m pixel with a Poisson distribution ---------------------- method 3
      # dbh$w25m <- dbh$wlid / 16
      # dbh$w25m <- rpois(nrow(dbh), dbh$w25m)

      # calculate weight for a 25*25m pixel depending on decimal of wlid ----------------------- method 4
      dbh$w25m <- dbh$wlid / 16
      dbh$decimal <- dbh$w25m - floor(dbh$w25m)
      dbh$random <- runif(nrow(dbh), min = 0, max = 1)
      dbh$randSmallerThanw25m <- dbh$random < dbh$decimal
      dbh[dbh$randSmallerThanw25m == TRUE, 'add'] <- 1
      dbh[dbh$randSmallerThanw25m == FALSE, 'add'] <- 0
      dbh$w25m <- floor(dbh$w25m) + dbh$add

      # remove trees with w25m = 0
      dbh <- dbh[dbh$w25m > 0,]

      # assign new dbh to trees while keeping BAtreelid
      dbh$dbhlid25m <-sqrt(40000/pi * (dbh$BAtreelid/16) / dbh$w25m)

      # return
      df <- dbh[, c('species_name', 'wlid', 'w25m', 'dbhlid25m', 'i', 'cellID25')]
      df$wlid <- df$wlid/16
      colnames(df) <- c('sp', 'wlid', 'n', 'dbh', 'i', 'cellID25')
      # if all trees are removed because of their small weights
      # an empty df must be created
      if(nrow(df) == 0){
        df <- data.frame(sp = NA, wlid = NA, n = NA, dbh = NA, i = i, cellID25 = cellID25)
      }

    } else{
      df <- data.frame(sp = NA, wlid = NA, n = NA, dbh = NA, i = i, cellID25 = cellID25)
    }

    return(df)

  }

  # parallel calculation on raster cells (one raster after the other)
  clustCalc <- function(rast, assignDendro, tree, NFIsp){
    # set cluster
    cl <- makeCluster(5)
    registerDoParallel(cl)
    results <- foreach(i = 1:nrow(rast[]), .combine = 'rbind', .packages = c('raster', 'rgdal')) %dopar% {assignDendro(cell = rast[i], i = i, tree, NFIsp)}
    stopCluster(cl)
    return(results)
  }

  # run calculation
  start <- Sys.time()
  results <- lapply(c(rast1, rast2, rast3, rast4), clustCalc, assignDendro, tree, NFIsp)
  end <- Sys.time()
  end - start

  # assemble results into one dataframe
  results1 <- as.data.frame(results[1])
  results2 <- as.data.frame(results[2])
  results3 <- as.data.frame(results[3])
  results4 <- as.data.frame(results[4])
  results2$i <- results2$i + max(results1$i)
  results3$i <- results3$i + max(results2$i)
  results4$i <- results4$i + max(results3$i)
  results <- rbind(results1, results2, results3, results4)

  # remove trees with sp = XXXX and n = dbh = NA
  # but only if it does not remove the cell
  # these 'strange' trees are created because they are deci or coni on cells
  # with dprop = 0 or 1
  # first count number of cells with composition = XXXX and with all trees dbh = NA
  df2 <- results %>% filter(!is.na(sp)) %>% group_by(cellID25) %>% summarise(N = sum(n, na.rm = TRUE))
  nrow(df2[df2$N == 0,]) # if equal 0, then no cells are removed only deci or coni trees in cells where dprop = 1 or 0
  # second delete trees with dbh = n = NA if it does not remove the cell
  if (nrow(df2[df2$N == 0,]) == 0){
    results[is.na(results$dbh), 'sp'] <- NA
  }

  # save
  saveRDS(results[, c('cellID25', 'sp', 'n', 'dbh', 'wlid', 'i')], file = paste0(tempPath, '/trees.rds'))

}
