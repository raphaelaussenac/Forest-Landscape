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
  require(ggplot2)

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


  ###############################################################
  # load tree data and vegetation type data
  ###############################################################

  # load TFV spatial data
  bd <- readOGR(dsn = './data/bauges/GEO', layer = 'BD_Foret_V2_PNRfilled_Foret_2014', encoding = 'UTF-8', use_iconv = TRUE)
  # load NFI tree data
  tree <- read.csv('./data/bauges/NFI/arbres_Bauges_2020_10_15.csv', sep = ';')
  # load correspondence between NFI plots and TFV types
  crpdTFV <- read.csv('./data/bauges/NFI/codeTFV_Bauges_2020_10_15.csv', sep = ';')


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
  # not specified --------------------------------------------> removed

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
  # create TFV raster
  ###############################################################

  # first assign TFV code to each forest cell
  bd$CODE_TFV <- factor(bd$CODE_TFV, levels = c('FF1-00-00', 'FF2G61-61', 'FF31', 'FF32', 'FF1-09-09'))
  bd$CODE_TFV <- as.numeric(bd$CODE_TFV)
  TFVraster <- rasterize(bd, dg, field = 'CODE_TFV')
  writeRaster(TFVraster, filename = paste0(tempPath, '/tfv.asc'), format = 'ascii', overwrite = TRUE)


  ###############################################################
  # calculate DBH and retrieve sp latin names
  ###############################################################

  # import latin names and create deciduous / coniferous categories
  tree <- spTransform(tree)

  # calculate proportion of deciduous BA
  tree$DBH <- tree$c13 / pi
  # save tree file for evaluation
  saveRDS(tree, file = paste0(tempPath, '/treeTemp.rds'))

}
