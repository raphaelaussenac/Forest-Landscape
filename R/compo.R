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
  require(doFuture)
  plan(multisession)

  # load LIDAR rasters
  Dg <- raster(paste0(tempPath, '/dg.grd'))
  BA <- raster(paste0(tempPath, '/BA.grd'))
  Dprop <- raster(paste0(tempPath, '/Dprop.grd'))
  if(landscape == 'bauges'){
    TFVraster <- raster(paste0(tempPath, '/tfv.asc'))
    names(TFVraster) <- 'CODE_TFV'
  }

  # load cellID raster
  cellID25 <- raster(paste0(landPath, '/cellID25.asc'))

  # load tree data and vegetation type data
  tree <- readRDS(paste0(tempPath, '/treeTemp.rds'))

  ###############################################################
  # Retrieve Dg, BA, Dprop and TFV for all plots
  ###############################################################

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
  # total range defined on NFI + LiDAR values
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
  fig <- plot_ly(x = NFI$Dg01, y = NFI$Dprop, z = NFI$BA01,
                 text = NFI$idp, textposition = 'middle right', # text
                 textfont = list(color = '#000000', size = 16), # text
                 type='scatter3d', mode = 'markers+text', color = NFI$CODE_TFV)
  axx <- list(title = 'quadratic diameter (cm)')
  axy <- list(title = 'proportion of Deciduous basal area')
  axz <- list(title = 'basal area (m^2/ha)')
  fig <- fig %>% layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  fig

  # convert NFI TFV CODE in numeric values
  if (landscape == 'bauges'){
      NFI$CODE_TFV <- factor(NFI$CODE_TFV, levels = c('FF1-00-00', 'FF2G61-61', 'FF31', 'FF32', 'FF1-09-09'))
  }
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
  # create CODE_TFV raster if not provided
  ###############################################################

  if(landscape != 'bauges'){
    TFVraster <- BA
    TFVraster[!is.na(TFVraster)] <- unique(NFI$CODE_TFV)
    names(TFVraster) <- 'CODE_TFV'
  }


  ###############################################################
  # Split raster into horizontal strips (speed-up calculation)
  ###############################################################

  # Split data into two separate stack rasters
  compoRaster <- stack(TFVraster, Dprop, Dg01, BA01, compo, cellID25)

  # gest compoRaster coordinates
  xmin <- xmin(compoRaster)
  xmax <- xmax(compoRaster)
  ymin <- ymin(compoRaster)
  ymax <- ymax(compoRaster)

  # nb of horizontal strips
  nbS <- 20

  # define ymin and ymax of horizontal strips
  ystep <- seq(from = ymin, to = ymax, by = (ymax - ymin) / nbS)

  strips <- list()
  # create nbS rasters
  for (i in 1:nbS){
    stripX <- crop(compoRaster, c(xmin, xmax, ystep[i], ystep[i+1]))
    strips <- c(strips, stripX)
  }


  ###############################################################
  # claculate distance between forest cells values and NFI values
  # and assign composition of nearest NFI plot
  ###############################################################

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
    cl <- makeCluster(6)
    registerDoParallel(cl)
    # results <- foreach(i = 1:nrow(rast[]), .combine = 'rbind', .packages = c('raster', 'rgdal')) %dopar% {assignCompo(cell = rast[i], i = i, NFI)}
    results <- foreach(i = 1:nrow(rast[]), .combine = 'rbind', .options.future = list(packages = c('raster', 'rgdal'))) %dofuture% {assignCompo(cell = rast[i], i = i, NFI)}
    stopCluster(cl)
    # transfer results into raster stack
    results <- data.frame(results)
    colnames(results) <- c('i', 'id')
    results <- results %>% arrange(i)
    rast$compo <- results[,2]
    # plot(rast$compo)
    return(rast)
  }

  # run calculation
  start <- Sys.time()
  rasts <- lapply(strips, clustCalc, assignCompo, NFI)
  end <- Sys.time()
  end - start

  ################################################################################
  # after splitting landscapes into strips (in dendro.R)
  # check whether all cells are present only once
  # check whether there is no missing cell
  ################################################################################

  # ensure cells are only present once (only in one raster strip)
  # get cellID25 of all strips
  getCellID <- function(raster){
    df <- as.data.frame(raster)
    cells <- unique(df$cellID25)
    return(cells)
  }
  cells <- lapply(rasts, getCellID)
  # are there duplicated cells
  cellList <- do.call(c, cells)
  dup <- sum(duplicated(cellList))

  # print error message if cells are duplicated
  if(dup > 0){
    stop("When splitting rasters into strips (in dendro.R) some cells are duplicated. Try to split ratser in a different number of strips")
  }
  # print error message if cells are missing
  if(length(cellList) != (nrow(cellID25) * ncol(cellID25))){
    stop("When splitting rasters into strips (in dendro.R) some cells are lost. Try to split ratser in a different number of strips")
  }

  ################################################################################
  # format and save compo raster
  ################################################################################

  # keep only compo in strips
  keepCompo <- function(raster){
    raster <- raster$compo
    return(raster)
  }
  compo <- lapply(rasts, keepCompo)

  # merge strips
  compo <- do.call(merge, compo)
  names(compo) <- 'compo'

  # save
  writeRaster(compo, paste0(tempPath, '/compoID.asc'), overwrite = TRUE)

}
