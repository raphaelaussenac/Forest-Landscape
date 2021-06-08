evalDendro <- function(){

  ###############################################################
  # initialisation
  ###############################################################

  # load packages
  require(dplyr)
  require(ggplot2)
  require(raster)

  # load tree data
  results <- read.csv(paste0(tempPath, '/trees75ForEval.csv'))

  ################################################################################
  # evaluate the effect of rounding on the number of trees
  ################################################################################

  # distribution of differences between rounded and non-rounded number of trees
  # at the cell level
  Ncell <- results %>% group_by(cellID25) %>% summarise(wlid = sum(wlid), n = sum(n)) %>%
                       mutate(Ndiff = wlid - n)
  pl1 <- ggplot()+
  geom_histogram(data = Ncell, aes(x = Ndiff), bins = 50, col = 'black', fill = 'white', alpha = 0.5) +
  theme_bw() +
  xlab('N diff at the cell level (negative values = overestimation)')
  pl1
  ggsave(file = paste0(evalPath, '/NdiffCell.pdf'), plot = pl1, width = 10, height = 10)

  # distribution of differences between rounded and non-rounded number of trees per species
  # at the cell level depending on dbh classes
  sp <- results %>% group_by(cellID25, sp) %>% summarise(SpMeanDBH = sum(n * dbh) / sum(n), SpNDiff = sum(wlid) - sum(n))
  sp$meanDBHclass <- cut(sp$SpMeanDBH, 0:120)
  pl2 <- ggplot() +
  geom_boxplot(data = sp, aes(x = meanDBHclass, y = SpNDiff), alpha = 0.5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab('species mean dbh in cells') +
  ylab('species N diff at the cell level (negative values = overestimation)')
  pl2
  ggsave(file = paste0(evalPath, '/NdiffDBH.pdf'), plot = pl2, width = 30, height = 20)

  # distribution of differences between rounded and non-rounded number of trees per species
  # at the cell level depending on dbh classes separately for each species
  pl3 <- ggplot() +
  geom_boxplot(data = sp, aes(x = meanDBHclass, y = SpNDiff), alpha = 0.5) +
  facet_wrap(~ sp) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(colour = 'black'),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text=element_text(size = 3)) +
  xlab('species mean dbh in cells') +
  ylab('species N diff at the cell level (negative values = overestimation)')
  pl3
  ggsave(file = paste0(evalPath, '/NdiffDBHSp.pdf'), plot = pl3, width = 30, height = 20)

  # total over/under estimation of trees at the landscape scale
  # positive values = underestimation
  sum(Ncell$Ndiff, na.rm = TRUE)
  # relative to the total number of trees (in %)
  sum(Ncell$Ndiff, na.rm = TRUE) * 100 /sum(results$n, na.rm = TRUE)

  ################################################################################
  # check whether BAlidar = sum of BA of trees at each cell
  ################################################################################

  BA <- raster(paste0(tempPath, '/BA.asc'))
  BA[BA > 1000] <- 1000 # correct the wrong min max values

  BAinit <- results %>% group_by(cellID25) %>% summarise(BAtot = sum((pi * (dbh/200)^2) * n * 16)) %>% arrange(cellID25)
  BA$BAinit <- BAinit$BAtot
  BA$BAdiff25m <- (BA$BA / 16)  - (BA$BAinit /16)
  BA$BAreldiff25m <- 100 - ( (BA$BAinit /16) * 100 / (BA$BA / 16) )

  # difference should be as close to zero as possible
  pdf(file = paste0(evalPath, '/BAdiff25m.pdf'))
  hist(BA$BAdiff25m, breaks = 100)
  dev.off()
  plot(BA$BAdiff25m)
  hist(BA$BAreldiff25m, breaks = 100)
  plot(BA$BAreldiff25m)

  ################################################################################
  # check whether Dglidar = mean Dg of trees at each cell
  ################################################################################

  Dg <- raster(paste0(tempPath, '/dg.asc'))
  Dg[Dg > 1000] <- 1000 # correct the wrong min max values

  Dginit <- results %>% group_by(cellID25) %>% summarise(Dgtot = sqrt(sum(dbh^2 * n)/sum(n))) %>% arrange(cellID25)
  Dg$Dginit <- Dginit$Dgtot
  Dg$Dgdiff <- Dg$dg - Dg$Dginit
  Dg$Dgreldiff <- 100 - (Dg$Dginit * 100 /Dg$dg)

  # difference should be as close to zero as possible but it necessarily
  # fluctuates because we had to change the trees dbh to reach the BA
  # prescribed by the LIDAR
  pdf(file = paste0(evalPath, '/Dgdiff.pdf'))
  hist(Dg$Dgdiff, breaks = 100)
  dev.off()
  plot(Dg$Dgdiff)
  hist(Dg$Dgreldiff, breaks = 100)
  plot(Dg$Dgreldiff)
  # writeRaster(Dg$Dgreldiff, filename = './data/Init/rasterVerif.asc', format = 'ascii', overwrite = TRUE)
  # Dg[802041]
  # results[results$cellID == 802041,]
  # saveResults[saveResults$cellID == 802041,]
  #
  # # check forest surface
  # sum(!is.na(rasterStack$dg[])) * 625 / 10000
  # sum(!is.na(rasterStack$compoID[])) * 625 /10000
  # sum(!is.na(rasterStack$BAinit[])) * 625 / 10000
  # sum(!is.na(rasterStack$Dginit[])) * 625 / 10000


}
