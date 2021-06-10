evalCompo <- function(){

  ################################################################################
  # initialisation
  ################################################################################

  # load packages
  require(dplyr)
  require(raster)
  require(rgdal)
  require(ggplot2)

  # load treeTemp
  # tree <- read.csv(paste0(tempPath, '/treeTemp.csv'))
  tree <- readRDS(paste0(tempPath, '/treeTemp.rds'))

  # load compoID
  compo <- raster(paste0(tempPath, '/compoID.asc'))


  ################################################################################
  # count number of times NFI plots are picked to define cells compo
  ################################################################################

  compodf <- as.data.frame(compo)
  compodf$compoID <- as.factor(compodf$compoID)
  compodf <- compodf %>% filter(!is.na(compoID)) %>%
                     group_by(compoID) %>%
                     summarise(count = n()) %>%
                     arrange(-count)
  compodf$compoID <- factor(compodf$compoID, levels = compodf$compoID)

  # plot
  plpl <- ggplot(data = compodf, aes(x = compoID, y = count)) +
  geom_bar(stat = 'identity') +
  xlab('plot ID') +
  ggtitle('number of times each plot is picked to define cells composition') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  plpl
  ggsave(file = paste0(evalPath, '/plotDraw.pdf'), plot = plpl, width = 25, height = 10)


  ################################################################################
  # check composition - convert plot id map into species composition map
  ################################################################################

  # set threshold to identify mainSp
  # thresh = 0.8 means you will get the species making up for > 80% of the stand BA
  # the stand BA
  thresh = 0.70
  # retrieve main species
  mainSp <- tree %>% group_by(idp, species_name) %>%
                            summarise(BA = sum((pi * (DBH/200)^2) * w)) %>%
                            group_by(idp) %>% arrange(idp, -BA) %>%
                            mutate(BAtot = sum(BA), BAprop = BA/BAtot, cumulProp = cumsum(BAprop), HigherThanthresh = case_when(
                              cumulProp >= thresh ~ cumulProp), minCompo = min(HigherThanthresh, na.rm = TRUE)) %>%
                            filter(cumulProp <= minCompo) %>% summarise(sp = paste(species_name, collapse=' - '))
  #
  # convert raster into polygon
  compoPoly <- rasterToPolygons(compo$compoID, n = 4, na.rm = TRUE, digits=12, dissolve = TRUE)
  # transfer main species values into the polygon
  compoPoly <- merge(compoPoly, mainSp, by.x = 'compoID', by.y = 'idp')
  writeOGR(compoPoly, evalPath, 'compoPoly', driver = 'ESRI Shapefile', overwrite = TRUE)

  # calculate surface for each stand type (based on their main sp)
  compoPoly$area <- area(compoPoly)
  surf <- as.data.frame(compoPoly) %>% group_by(sp) %>% summarise(Area = sum(area)/10000) %>% arrange(-Area)
  surf$sp <- factor(surf$sp, levels = as.character(surf$sp))
  # plot results
  pl1 <- ggplot(data = surf) +
  geom_bar(aes(x = sp, y = Area), stat = 'identity') +
  theme_light() +
  xlab('main species') +
  ylab('surface (ha)') +
  # scale_y_continuous('surface (ha)', sec.axis = sec_axis(~ . * 1.20, name = 'number of NFI plots')) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        plot.title = element_text(hjust = 0.5))
  pl1
  ggsave(file = paste0(evalPath, '/compoSurface.pdf'), plot = pl1, width = 20, height = 10)

}
