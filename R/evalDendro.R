evalDendro <- function(){

  ###############################################################
  # initialisation
  ###############################################################

  # load packages
  require(dplyr)
  require(ggplot2)
  require(raster)

  # load tree data
  results <- readRDS(paste0(tempPath, '/trees75.rds'))

  # load lidar data
  Dg <- raster(paste0(tempPath, '/dg.grd'))
  BA <- raster(paste0(tempPath, '/BA.grd'))
  cellID25 <- raster(paste0(landPath, '/cellID25.asc'))
  rasterStack <- stack(Dg, BA, cellID25)
  df <- as.data.frame(rasterStack)


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
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'red') +
  xlab('N diff at the cell level (negative values = overestimation)')
  pl1
  ggsave(file = paste0(evalPath, '/NdiffCell.pdf'), plot = pl1, width = 10, height = 10)

  # distribution of differences between rounded and non-rounded number of trees
  # at the cell level depending on dbh classes
  Ncelldbh <- results %>% group_by(cellID25) %>% summarise(meanDBH = sum(n * dbh) / sum(n), NDiff = sum(wlid) - sum(n))
  Ncelldbh$meanDBHclass <- cut(Ncelldbh$meanDBH, 0:150)
  pl2 <- ggplot() +
  geom_boxplot(data = Ncelldbh, aes(x = meanDBHclass, y = NDiff), alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'red') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab('mean dbh in cells') +
  ylab('N diff at the cell level (negative values = overestimation)')
  pl2
  ggsave(file = paste0(evalPath, '/NdiffDBH.pdf'), plot = pl2, width = 30, height = 20)

  # distribution of differences between rounded and non-rounded number of trees per species
  # at the cell level depending on dbh classes separately for each species
  sp <- results %>% group_by(cellID25, sp) %>% summarise(SpMeanDBH = sum(n * dbh) / sum(n), SpNDiff = sum(wlid) - sum(n))
  sp$meanDBHclass <- cut(sp$SpMeanDBH, 0:150)
  pl3 <- ggplot() +
  geom_boxplot(data = sp, aes(x = meanDBHclass, y = SpNDiff), alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'red') +
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
  # quantiles of Ndiff
  quantile(Ncell$Ndiff)

  ################################################################################
  # check whether BAlidar = sum of BA of trees at each cell
  ################################################################################

  BAinit <- results %>% group_by(cellID25) %>% summarise(BAtot = sum((pi * (dbh/200)^2) * n)) %>% arrange(cellID25)
  df <- merge(df, BAinit, by = 'cellID25', all.x = TRUE)
  df$BAdiff25m <- (df$BA / 16) - df$BAtot
  df$BAreldiff25m <- 100 - ( df$BAtot * 100 / (df$BA / 16) )

  # difference should be as close to zero as possible
  pdf(file = paste0(evalPath, '/BAdiff25m.pdf'))
  hist(df$BAdiff25m, breaks = 100)
  dev.off()
  hist(df$BAreldiff25m, breaks = 100)
  rasterStack$BAdiff25m <- df$BAdiff25m
  rasterStack$BAreldiff25m <- df$BAreldiff25m
  plot(rasterStack$BAdiff25m)
  plot(rasterStack$BAreldiff25m)

  ################################################################################
  # check whether Dglidar = mean Dg of trees at each cell
  ################################################################################

  Dginit <- results %>% group_by(cellID25) %>% summarise(Dgtot = sqrt(sum(dbh^2 * n)/sum(n))) %>% arrange(cellID25)
  df <- merge(df, Dginit, by = 'cellID25', all.x = TRUE)
  df$Dgdiff <- df$Dg - df$Dgtot
  df$Dgreldiff <- 100 - (df$Dgtot * 100 /df$Dg)

  # difference should be as close to zero as possible but it necessarily
  # fluctuates because we had to change the trees dbh to reach the BA
  # prescribed by the LIDAR
  pdf(file = paste0(evalPath, '/Dgdiff.pdf'))
  hist(df$Dgdiff, breaks = 100)
  dev.off()
  hist(df$Dgreldiff, breaks = 100)
  rasterStack$Dgdiff <- df$Dgdiff
  rasterStack$Dgreldiff <- df$Dgreldiff
  plot(rasterStack$Dgdiff)
  plot(rasterStack$Dgreldiff)
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
