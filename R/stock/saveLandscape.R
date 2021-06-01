saveLandscape <- function(cellID25, results){

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
  writeRaster(cellID100, filename = './initialLandscape/cellID100.asc', format = 'ascii', overwrite = TRUE)
  # change resolution back to cellID25 resolution in order to crop and then
  # stack together
  cellID100 <- disaggregate(cellID100, fact = resFact)
  # reduce extent of cellID100 to minimum common extent with cellID25
  cellID100 <- crop(cellID100, extent(cellID25))
  # stack rasters
  cellID <- stack(cellID25, cellID100)
  # cellID correspondence dataframe
  cellIDdf <- as.data.frame(values(cellID))
  # merge with other results
  results <- merge(results, cellIDdf, by = 'cellID25')

  ################################################################################
  # save
  ################################################################################

  # add cellID100 to environmental data
  envdf <- read.csv('./data/init/envVariablesTemp.csv')
  envdf <- merge(envdf, cellIDdf, by = 'cellID25', all.y = F)
  # nb of forest cells per ha
  envdf <- envdf %>% group_by(cellID100) %>% mutate(forestCellsPerHa = sum(forest))
  # reduce table size in memory
  envdf$cellID100 <- as.integer(envdf$cellID100)
  # sort colnames
  colOrd <- c('cellID25','cellID100','park','forest','elev','slope','aspect','swhc','pH','GRECO','SIQpet','SIFsyl','SIAalb','SIPabi','forestCellsPerHa')
  # save
  write.csv(envdf[, colOrd], file = './initialLandscape/envVariables.csv', row.names = FALSE)

  # remove cells with no trees
  results <- results[!is.na(results$dbh),]
  # reduce table size in memory
  results$cellID25 <- as.integer(results$cellID25)
  results$cellID100 <- as.integer(results$cellID100)
  results$n <- as.integer(results$n)
  results$dbh <- round(results$dbh, 2)

  # save
  write.csv(results[, c('cellID25', 'cellID100', 'sp', 'n', 'dbh')], file = './initialLandscape/trees75.csv', row.names = FALSE)

  return(results)

}
