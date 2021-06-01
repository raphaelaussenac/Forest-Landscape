compo <- function(landscape){

  ###############################################################
  # initialisation
  ###############################################################

  # load packages
  require(rgdal)
  require(raster)
  require(dplyr)
  require(ggplot2)
  require(plotly)
  require(doParallel)

  # load LIDAR rasters
  Dg <- raster(paste0(tempPath, '/dg.asc'))
  Dg[Dg > 1000] <- 1000 # correct the wrong min max values
  BA <- raster(paste0(tempPath, '/BA.asc'))
  BA[BA > 1000] <- 1000 # correct the wrong min max values
  Dprop <- raster(paste0(tempPath, '/Dprop.asc'))
  Dprop <- Dprop / 100

  # load tree data and vegetation type data
  if(landscape == 'bauges'){
    # load TFV spatial data
    bd <- readOGR(dsn = './data/bauges/GEO', layer = 'BD_Foret_V2_PNRfilled_Foret_2014', encoding = 'UTF-8', use_iconv = TRUE)
    # load NFI tree data
    tree <- read.csv('./data/bauges/NFI/arbres_Bauges_2020_10_15.csv', sep = ';')
    # load correspondence between NFI plots and TFV types
    crpdTFV <- read.csv('./data/bauges/NFI/codeTFV_Bauges_2020_10_15.csv', sep = ';')
  }


  ###############################################################
  # group TFV types together
  ###############################################################

  # first assign TFV to each NFI plot using crpdTFV
  tree <- merge(tree, crpdTFV[, c('idp', 'tfv')], by = 'idp', all.x = TRUE)
  tree$tfv <- as.factor(tree$tfv)
  colnames(tree)[colnames(tree) == 'tfv'] <- 'CODE_TFV'

  # Surface area of each TFV type and number of associated NFI plots
  TFVcountAndSurface <- function(bd, tree){
    # calculate area of each TFV type
    bd$area <- area(bd)
    surfArea <- data.frame(bd) %>% group_by(CODE_TFV) %>%
                                  summarize(surface = sum(area)/10000) %>%
                                  arrange(-surface)
    # count number of NFI in each TFV type
    TFVplotCount <- tree %>% group_by(CODE_TFV) %>% summarise(N = length(unique(idp)))
    # merge with TFV surface
    surfArea <- merge(surfArea, TFVplotCount, by = 'CODE_TFV', all = TRUE)
    surfArea <- surfArea %>% arrange(-surface)
    surfArea[is.na(surfArea$surface), 'surface'] <- 0
    surfArea$CODE_TFV <- factor(surfArea$CODE_TFV, levels = as.character(surfArea$CODE_TFV))

    # first five TFV types represent 96% of the total forest cover
    prop <- sum(surfArea[1:5, 'surface']) * 100 / sum(surfArea$surface)

    # plot
    pl1 <- ggplot(data = surfArea) +
    geom_bar(aes(x = CODE_TFV, y = surface), stat = 'identity') +
    geom_bar(data = surfArea[1:5,], aes(x = CODE_TFV, y = surface), stat = 'identity', col = 'black', fill = 'black') +
    geom_bar(aes(x = CODE_TFV, y = N * 20), stat = 'identity', col = 'orange', fill = 'orange', alpha = 0.5) +
    geom_text(aes(x = CODE_TFV, y = N * 20, label = N), vjust = -0.3, size = 3.5, col = 'orange') +
    geom_text(data = surfArea[3,], aes(x = CODE_TFV, y = surface, label =  paste(round(prop, 2), '% of forest cover')), vjust = -1, hjust = 0, size = 3.5, col = 'black') +
    # annotate('text', x = 3, y = 10000, label = paste(round(prop, 2), '% of forest cover')) +
    theme_light() +
    xlab('TFV code') +
    ylab('surface (ha)') +
    # scale_y_continuous('surface (ha)', sec.axis = sec_axis(~ . * 1.20, name = 'number of NFI plots')) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          plot.title = element_text(hjust = 0.5))

    return(pl1)

  }

  pl1 <- TFVcountAndSurface(bd, tree)
  pl1 + ggtitle('TFV surface and number of NFI plots (before grouping TFV types together)')

  # group TFV types together
  # main types
  A <- 'FF1-00-00' #: Forêt fermée à mélange de feuillus
  B <- 'FF2G61-61' #: Forêt fermée de sapin ou épicéa
  C <- 'FF31' #: Forêt fermée à mélange de feuillus prépondérants et conifères
  D <- 'FF32' #: Forêt fermée à mélange de conifères prépondérants et feuillus
  E <- 'FF1-09-09' #: Forêt fermée de hêtre pur

  # under-represented types
  # FO1: Forêt ouverte de feuillus purs ----------------------> A
  # FF1-00: Forêt fermée de feuillus purs en îlots -----------> A
  # FO2: Forêt ouverte de conifères purs ---------------------> B
  # FO3: Forêt ouverte à mélange de feuillus et conifères ----> C
  # FF1-49-49: Forêt fermée d’un autre feuillu pur -----------> A
  # FF2-00-00: Forêt fermée à mélange de conifères -----------> B
  # FF1G01-01: Forêt fermée de chênes décidus purs -----------> A
  # FF2G53-53: Forêt fermée de pin laricio ou pin noir pur ---> B
  # FF2-90-90: Forêt fermée à mélange d’autres conifères -----> B
  # FF2-64-64: Forêt fermée de douglas pur -------------------> B
  # FF2-00: Forêt fermée de conifères purs en îlots ----------> B
  # FF1-10-10: Forêt fermée de châtaignier pur ---------------> A
  # FO0: Forêt ouverte sans couvert arboré -------------------> remove
  # FF2-52-52: Forêt fermée de pin sylvestre pur -------------> B
  # FF0: Forêt fermée sans couvert arboré --------------------> remove
  # LA4: Formation herbacée ----------------------------------> remove
  # LA6: Lande -----------------------------------------------> remove
  # not specified --------------------------------------------> A

  groupTfv <- function(df){
    df <- df[!(df$CODE_TFV %in% c('FO0', 'FF0', 'LA4', 'LA6')) & !is.na(df$CODE_TFV),]
    df[df$CODE_TFV == 'FO1', 'CODE_TFV'] <- A
    df[df$CODE_TFV == 'FF1-00', 'CODE_TFV'] <- A
    df[df$CODE_TFV == 'FO2', 'CODE_TFV'] <- B
    df[df$CODE_TFV == 'FO3', 'CODE_TFV'] <- C
    df[df$CODE_TFV == 'FF1-49-49', 'CODE_TFV'] <- A
    df[df$CODE_TFV == 'FF2-00-00', 'CODE_TFV'] <- B
    df[df$CODE_TFV == 'FF1G01-01', 'CODE_TFV'] <- A
    df[df$CODE_TFV == 'FF2G53-53', 'CODE_TFV'] <- B
    df[df$CODE_TFV == 'FF2-90-90', 'CODE_TFV'] <- B
    df[df$CODE_TFV == 'FF2-64-64', 'CODE_TFV'] <- B
    df[df$CODE_TFV == 'FF2-00', 'CODE_TFV'] <- B
    df[df$CODE_TFV == 'FF1-10-10', 'CODE_TFV'] <- A
    df[df$CODE_TFV == 'FF2-52-52', 'CODE_TFV'] <- B
    df[df$CODE_TFV == '', 'CODE_TFV'] <- A
    return(df)
  }

  # assign new TFV to spatial polygons
  bd$CODE_TFV <- as.character(bd$CODE_TFV)
  bd <- groupTfv(bd)

  # assign new TFV to NFI data
  tree$CODE_TFV <- as.character(tree$CODE_TFV)
  tree <- groupTfv(tree)

  pl1 <- TFVcountAndSurface(bd, tree)
  pl1 + ggtitle('TFV surface and number of NFI plots (after grouping TFV types together)')


  ###############################################################
  # Retrieve Dg, BA, Dprop and TFV for all forest cells
  ###############################################################

  # first assign TFV code to each forest cell
  bd$CODE_TFV <- factor(bd$CODE_TFV, levels = c('FF1-00-00', 'FF2G61-61', 'FF31', 'FF32', 'FF1-09-09'))
  bd$CODE_TFV <- as.numeric(bd$CODE_TFV)
  TFVraster <- rasterize(bd, Dg, field = 'CODE_TFV')
  names(TFVraster) <- 'CODE_TFV'

  ###############################################################
  # Retrieve Dg, BA, Dprop and TFV for all NFI plots
  ###############################################################

  # import latin names and create deciduous / coniferous categories
  tree <- spTransform(tree)

  # calculate proportion of deciduous BA
  tree$DBH <- tree$c13 / pi
  # save tree file for evaluation
  write.csv(tree, file = paste0(tempPath, '/treeTemp.csv'), row.names = FALSE)
  NFIDprop <- tree %>% group_by(idp, spType) %>%
                            summarise(BAdc = sum((pi * (DBH/200)^2) * w)) %>%
                            group_by(idp) %>% mutate(BA = sum(BAdc), Dprop = BAdc/BA) %>%
                            filter(spType == 'D') %>% select(-spType, -BAdc, -BA)
  # calculate Dg and BA
  NFI <- tree %>% group_by(idp) %>%
                      summarise(Dg = sqrt(sum(DBH^2 * w)/sum(w)),
                                BA = sum((pi * (DBH/200)^2) * w),
                                CODE_TFV = unique(CODE_TFV))
  NFI <- merge(NFIDprop, NFI, by = 'idp', all = TRUE)
  NFI[is.na(NFI$Dprop), 'Dprop'] <- 0

  ###############################################################
  # normalise Dg and BA values to 0-1 range
  # range defined on NFI values
  ###############################################################

  # NFI plots
  minDg <- min(min(NFI$Dg), minValue(Dg))
  maxDg <- max(max(NFI$Dg), maxValue(Dg))
  NFI$Dg01 <- ( NFI$Dg - minDg ) / ( maxDg - minDg )
  minBA <- min(min(NFI$BA), minValue(BA))
  maxBA <- max(max(NFI$BA), maxValue(BA))
  NFI$BA01 <- ( NFI$BA - minBA ) / ( maxBA - minBA )

  # NFI[NFI$idp == '469390', 'CODE_TFV'] <- 'CASTA'

  # plot 3d
  fig <- plot_ly(x = NFI$Dg01, y = NFI$Dprop, z = NFI$BA01, type='scatter3d', mode = 'markers', color = NFI$CODE_TFV)
  axx <- list(title = 'quadratic diameter (cm)')
  axy <- list(title = 'proportion of Deciduous basal area')
  axz <- list(title = 'basal area (m^2/ha)')
  fig <- fig %>% layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  fig

  # convert NFI TFV CODE in numeric values
  NFI$CODE_TFV <- factor(NFI$CODE_TFV, levels = c('FF1-00-00', 'FF2G61-61', 'FF31', 'FF32', 'FF1-09-09'))
  NFI$CODE_TFV <- as.numeric(NFI$CODE_TFV)
  NFI <- NFI %>% select(idp, Dprop, Dg01, BA01, CODE_TFV)

  # forest cells
  Dg01 <- ( Dg - minDg ) / ( maxDg - minDg )
  names(Dg01) <- 'Dg01'
  BA01 <- ( BA - minBA ) / ( maxBA - minBA )
  names(BA01) <- 'BA01'

  ###############################################################
  # create empty raster of composition
  ###############################################################

  compo <- BA
  compo[!is.na(compo)] <- NA
  names(compo) <- 'compo'

  ###############################################################
  # claculate distance between forest cells values and NFI values
  # and assign composition of nearest NFI plot
  ###############################################################

  # Split data into two separate stack rasters
  compoRaster <- stack(TFVraster, Dprop, Dg01, BA01, compo)

  # rast1 <- crop(compoRaster, extent(compoRaster)/10)
  # rast2 <- crop(compoRaster, extent(compoRaster)/20)

  rast1 <- crop(compoRaster, c(extent(compoRaster)[1],
                                    extent(compoRaster)[1] + round( (extent(compoRaster)[2] - extent(compoRaster)[1]) / 2),
                                    extent(compoRaster)[3],
                                    extent(compoRaster)[4]))
  #
  rast2 <- crop(compoRaster, c(extent(compoRaster)[1] + round( (extent(compoRaster)[2] - extent(compoRaster)[1]) / 2),
                                    extent(compoRaster)[2],
                                    extent(compoRaster)[3],
                                    extent(compoRaster)[4]))
  #

  # function to assign composition
  assignCompo <- function(cell, i, NFI){
    # print(i)
    # some cells do not have any TFV value
    if(!is.na(cell[, 'CODE_TFV'])){
      # retrieve list of plots (and their Dg, BA, Dprop values) in the TFV
      # type associated to the cell
      plotlist <- NFI[NFI$CODE_TFV == cell[, 'CODE_TFV'], ]
      rownames(plotlist) <- plotlist$idp
      plotlist[, c('idp', 'CODE_TFV')] <- NULL
      # retrieve cell values = coordinates
      plotlist$DpropCell <- cell[, 'Dprop']
      plotlist$Dg01Cell <- cell[, 'Dg01']
      plotlist$BA01Cell <- cell[, 'BA01']
      # calculate distance
      plotlist$distance <- apply(plotlist, 1, function(x) dist(matrix(x, nrow = 2, byrow = TRUE)))
      # determine nearest NFI plot
      NearestPlot <- as.numeric(rownames(plotlist[plotlist$distance == min(plotlist$distance),]))
      NearestPlot <- NearestPlot[1]
    } else {
      NearestPlot <- NA
    }
    # assign nearest plot value (or NA) to cells
    return(c(i,NearestPlot))
  }

  # parallel calculation on raster cells (one raster after the other)
  clustCalc <- function(rast, assignCompo, NFI){
    # set cluster
    cl <- makeCluster(8)
    registerDoParallel(cl)
    results <- foreach(i = 1:nrow(rast[]), .combine = 'rbind', .packages = c('raster', 'rgdal')) %dopar% {assignCompo(cell = rast[i], i = i, NFI)}
    stopCluster(cl)
    # transfer results into raster stack
    results <- data.frame(results)
    colnames(results) <- c('i', 'id')
    results <- results %>% arrange(i)
    rast$compo <- results[,2]
    plot(rast)
    return(rast)
  }

  # run calculation
  start <- Sys.time()
  rasts <- lapply(c(rast1, rast2), clustCalc, assignCompo, NFI)
  end <- Sys.time()
  end - start

  # merge back rasters
  rast1 <- rasts[[1]]
  rast2 <- rasts[[2]]
  rast3 <- raster::merge(rast1, rast2, overlap = FALSE)
  names(rast3) <- names(compoRaster)
  # save
  writeRaster(rast3$compo, paste0(tempPath, '/compoID.asc'), overwrite = TRUE)

}
