managTable <- function(){

  ###############################################################
  # initialisation
  ###############################################################

  # load packages
  require(reldist)
  require(dplyr)
  require(rgdal)
  require(raster)
  require(ggplot2)
  library(stringr)

  # load virtual tree data
  tree <- read.csv(paste0(landPath, '/trees75.csv'))

  # load environmental data
  env <- read.csv(paste0(landPath, '/envVariables.csv'))

  # load protected areas
  rb <- readOGR(dsn = './data/bauges/GEO', layer = 'reserves_biologiques', encoding = 'UTF-8', use_iconv = TRUE)
  rn <- readOGR(dsn = './data/bauges/GEO', layer = 'reserves_naturelles', encoding = 'UTF-8', use_iconv = TRUE)

  # load park limits
  park <- readOGR(dsn = './data/bauges/GEO', layer = 'park', encoding = 'UTF-8', use_iconv = TRUE)

  # load ownership
  own <- readOGR(dsn = './data/bauges/GEO', layer = 'Foret_publique_dep73-74_2814_dissolve', encoding = 'UTF-8', use_iconv = TRUE)

  # load accessibility
  access <- raster('./data/bauges/GEO/PNRfilled_F.distance.tif')
  access[access > 100000] <- 0 # correct values>100000m

  # load cellID100 raster
  cellID100 <- raster(paste0(landPath, '/cellID100.asc'))

  # load species self-thinning boundary parameters
  rdiParam <- read.table('./data/valeursCoefficientsRdi.txt', header = T)

  # load list of deciduous and coniferous sp
  deciduousSp <- readRDS('./data/deciduousSp.rds')
  coniferousSp <- readRDS('./data/coniferousSp.rds')


  ###############################################################
  # define protected areas
  ###############################################################

  # set shp extent (not necessary since we use intersect just below
  # but fixes a bugue in package raster when rasterizing afterwards...)
  rb <- crop(rb, cellID100)
  rn <- crop(rn, cellID100)

  # keep only 'réserve biologique intégrale'
  # and exclude 'réserve biologique dirigée'
  rb <- rb[rb$code_r_enp == 'I',]

  # merge protected areas
  protect <- raster::union(rb, rn)
  protect$protection <- 1

  # select protected areas only in study area
  protect <- raster::intersect(protect, park)

  # convert into raster
  # use getCover to define proportion of each 100*100m cell covered by polygon
  testAntiBugue <- rasterize(protect, cellID100, getCover = TRUE)
  protect <- rasterize(protect, cellID100, getCover = TRUE)
  names(protect) <- 'protect'

  # stack with cellID100
  protect <- stack(cellID100, protect)

  # convert into dataframe
  protect <- as.data.frame(protect)

  # if protect >= 0.5 then most of the cell is covered by protected
  # area --> replace by 1. if protect < 0.5 --> replace by 0.
  df <- protect %>% mutate(protect = if_else(protect >= 0.5, 1, 0))


  ###############################################################
  # calculate Gini index on 1ha cells and define stand type
  # (evenaged / unevenaged)
  ###############################################################

  # if gini<0.45 --> evenaged stand, else --> unevenaged stand
  gini <- tree %>% group_by(cellID100) %>% mutate(ba = pi * dbh^2 / 4) %>%
                   summarise(gini = gini(x = ba, weights = n),
                             BA = sum((pi * (dbh/200)^2) * n),
                             Dg = sqrt(sum(dbh^2 * n)/sum(n)),
                             meanH = sum(h * n) /sum(n)) %>%
                   mutate(structure = if_else(gini < 0.45, 'even', 'uneven'))
  gini$structure <- as.factor(gini$structure)

  # add non-forest cells
  df <- merge(df, gini, by = 'cellID100', all.x = TRUE)

  ###############################################################
  # define composition type
  ###############################################################

  # set threshold to identify mainSp
  # thresh = 0.8 means you will get the species making up for >= 80% of
  # the stand BA
  thresh = 0.75
  # retrieve main species
  mainSp <- tree %>% group_by(cellID100, sp) %>%
                            summarise(BA = sum((pi * (dbh/200)^2) * n)) %>%
                            group_by(cellID100) %>% arrange(cellID100, -BA) %>%
                            mutate(BAtot = sum(BA), BAprop = BA/BAtot, cumulProp = cumsum(BAprop), HigherThanthresh = case_when(
                              cumulProp >= thresh ~ cumulProp), minCompo = min(HigherThanthresh, na.rm = TRUE)) %>%
                            filter(cumulProp <= minCompo) %>% summarise(mainSp = paste(sp, collapse=' - '))
  #

  # class main species into composition types:
  # identify deciduous and coniferous stands
  mainSp <- mainSp %>% group_by(cellID100) %>% mutate(D = sum(str_detect(mainSp, deciduousSp)),
                                                      C = sum(str_detect(mainSp, coniferousSp)),
                                                      compoType = if_else(D>0 & C>0, 'DC', if_else(D>0 & C==0, 'D', 'C')))
  #
  # beech
  mainSp[mainSp$mainSp == 'Fagus sylvatica' , 'compoType'] <- 'beech'
  # fir and/or spruce
  fs <- c('Abies alba', 'Picea abies', 'Abies alba - Picea abies', 'Picea abies - Abies alba')
  mainSp[mainSp$mainSp %in% fs , 'compoType'] <- 'fir and/or spruce'
  # mixed beech - fir and/or spruce
  m <- c('Fagus sylvatica - Abies alba','Abies alba - Fagus sylvatica',
         'Fagus sylvatica - Picea abies','Picea abies - Fagus sylvatica',
         'Fagus sylvatica - Abies alba - Picea abies',
         'Fagus sylvatica - Picea abies - Abies alba',
         'Abies alba - Picea abies - Fagus sylvatica',
         'Abies alba - Fagus sylvatica - Picea abies',
         'Picea abies - Fagus sylvatica - Abies alba',
         'Picea abies - Abies alba - Fagus sylvatica')
  mainSp[mainSp$mainSp %in% m , 'compoType'] <- 'beech with fir and/or spruce'

  # subdivide DC into DC and DC with fir and or spruce
  mainSp <- mainSp %>% mutate(DCfs = sum(str_detect(mainSp, c('Picea abies', 'Abies alba'))),
                              DCfsc = if_else(C > DCfs, 1, 0))
  mainSp[mainSp$compoType == 'DC' & mainSp$DCfs > 0 & mainSp$DCfsc == 0, 'compoType'] <- 'D with fir and/or spruce'
  mainSp[mainSp$compoType == 'DC' & mainSp$DCfs > 0 & mainSp$DCfsc == 1, 'compoType'] <- 'DC with fir and/or spruce'

  # subdivide C into C and C with fir and or spruce
  mainSp[mainSp$compoType == 'C' & mainSp$DCfs > 0, 'compoType'] <- 'C with fir and/or spruce'

  # add to df
  df <- merge(df, mainSp[, c('cellID100', 'compoType')], by = 'cellID100', all.x = TRUE)


  ###############################################################
  # ownership
  ###############################################################

  # retrieve ownership in study area
  own <- raster::intersect(own, park)

  # convert into raster
  # use getCover to define proportion of each 100*100m cell covered by polygon
  own <- rasterize(own, cellID100, getCover = TRUE)
  names(own) <- 'public'

  # stack with cellID100
  own <- stack(cellID100, own)

  # convert into dataframe
  own <- as.data.frame(own)

  # if public >= 0.5 then most of the cell is public --> replace by public.
  # if public < 0.5 --> replace by private.
  own <- own %>% mutate(owner = if_else(public >= 0.5, 'public', 'private'))

  # add to df
  df <- merge(df, own[, c('cellID100', 'owner')], by = 'cellID100')


  ###############################################################
  # access
  ###############################################################

  # aggregate access value from 5*5m to 100*100m
  access <- aggregate(access, fact = 20, fun = 'mean', na.rm = TRUE)

  # stack with cellID100
  access <- stack(cellID100, access)

  # convert into dataframe
  access <- as.data.frame(access)
  names(access) <- c('cellID100', 'dist')

  # if access > 2km --> not accessible
  access <- access %>% mutate(access = if_else(dist <= 2000, 1, 0)) %>%
                       dplyr::select(-dist)
  access[is.na(access$access), 'access'] <- 0
  # access = 1 --> accessible
  # access = 0 --> inaccessible

  # add to df
  df <- merge(df, access, by = 'cellID100')


  ###############################################################
  # number of forest 25*25m cells in 100*100m cells
  ###############################################################

  # nb of forest cells per ha
  forCel <- tree %>% group_by(cellID100) %>% summarise(forestCellsPerHa = length(unique(cellID25)))

  # add to df
  df <- merge(df, forCel, by = 'cellID100', all.x = TRUE)

  ###############################################################
  # calculate rdi
  ###############################################################

  # convert sp code into latinName
  rdiParam <- rdiParam %>% rename(espar = espece)
  rdiParam <- spTransform(rdiParam) %>% rename(sp = species_name)

  # number of tree of each sp on each cell
  rdi <- tree %>% group_by(cellID100, sp) %>% summarise(nsp = sum(n),
                                                        Dgsp = sqrt(sum(dbh^2 * n)/sum(n)))
  rdi <- merge(rdi, rdiParam, by = 'sp')

  # retrieve number of 25*25m forest cells in 100*100m cells
  rdi <- merge(rdi, df[, c('cellID100', 'forestCellsPerHa')], by = 'cellID100', all.x = TRUE)

  # calculate species Nmax
  rdi$Nspmax <- exp( rdi$rqIntercept + rdi$rqSlope * log(rdi$Dgsp) ) * (rdi$forestCellsPerHa / 16)

  # Calculate sp partial rdi
  rdi$rdip <- rdi$nsp / rdi$Nspmax

  # claculate total rdi
  rdit <- rdi %>% group_by(cellID100) %>% summarise(rdi = sum(rdip))

  # add to df
  df <- merge(df, rdit, by = 'cellID100', all.x = TRUE)

  # define density class for uneven-aged stands
  qtuneven <- quantile(df[df$structure == 'uneven' & df$access == 1, 'rdi'], na.rm = T, c(0.33, 0.66))
  hist(df[df$structure == 'uneven' & df$access == 1, 'rdi'], breaks = 100)
  abline(v = qtuneven, col = 'red', lty = 5, lwd = 2)
  df[df$structure == 'uneven' & !is.na(df$structure), 'density'] <- 'medium'
  df[df$structure == 'uneven' & !is.na(df$structure) & df$rdi > qtuneven[2], 'density'] <- 'high'
  df[df$structure == 'uneven' & !is.na(df$structure) & df$rdi < qtuneven[1], 'density'] <- 'low'

  # define density class for even-aged stands
  qteven <- quantile(df[df$structure == 'even' & df$access == 1, 'rdi'], na.rm = T, 0.5)
  hist(df[df$structure == 'even' & df$access == 1, 'rdi'], breaks = 100)
  abline(v = qteven, col = 'red', lty = 5, lwd = 2)
  df[df$structure == 'even' & !is.na(df$structure), 'density'] <- 'low'
  df[df$structure == 'even' & !is.na(df$structure) & df$rdi > qteven, 'density'] <- 'high'


  ###############################################################
  # define stand type (compoType + even / uneven / protect / ...)
  ###############################################################

  df$stand <- NA
  # no management (inaccessible + protected)
  df[(df$access == 0 | df$protect == 1) & !is.na(df$compoType), 'stand'] <- paste(df[(df$access == 0 | df$protect == 1) & !is.na(df$compoType), 'compoType'],
                                                              df[(df$access == 0 | df$protect == 1) & !is.na(df$compoType), 'owner'],
                                                              'noMan', sep = '-')
  # management
  df[(df$access == 1 & df$protect == 0) & !is.na(df$compoType), 'stand'] <- paste(df[(df$access == 1 & df$protect == 0) & !is.na(df$compoType), 'compoType'],
                                                              df[(df$access == 1 & df$protect == 0) & !is.na(df$compoType), 'owner'],
                                                              df[(df$access == 1 & df$protect == 0) & !is.na(df$compoType), 'structure'],
                                                              df[(df$access == 1 & df$protect == 0) & !is.na(df$compoType), 'density'], sep = '-')
  #
  ######################################################################################
  # assign 'final cut' management to XX% of some stand types
  # in priority to denser even-aged stands

  # fir and or spruce - private: 25%
  df1 <- df[df$owner == 'private' & df$compoType == 'fir and/or spruce' & !is.na(df$compoType), ] # 3932 ha en tout
  prop1 <- round(nrow(df1) * 0.25) # dont 983 devraient être 'final cut'
  df1 <- df1 %>% filter(structure == 'even', protect == 0, access == 1) %>% arrange(-rdi) # mais que 269 ha éligible (773 avec gini = 0.5)
  # cellSelect1 <- df1$cellID100[1:prop1]
  # df[df$cellID100 %in% cellSelect1, 'stand'] <- 'fir and/or spruce-private-final cut'

  # beech: 22%
  df2 <- df[df$compoType == 'beech' & !is.na(df$compoType), ] # 611
  prop2 <- round(nrow(df2) * 0.22) # 134
  df2 <- df2 %>% filter(structure == 'even', protect == 0, access == 1) %>% arrange(-rdi) # 144 (234 avec gini = 0.5)

  # other deciduous: 43%
  df3 <- df[df$compoType %in% c('D') & !is.na(df$compoType), ] # 10090
  prop3 <- round(nrow(df3) * 0.43) # 4339
  df3 <- df3 %>% filter(structure == 'even', protect == 0, access == 1) %>% arrange(-rdi) # 2224 (4343 avec gini = 0.5)

  # mixed - private: 4%
  df4 <- df[df$compoType %in% c('D with fir and/or spruce',
                                'beech with fir and/or spruce',
                                'DC with fir and/or spruce',
                                'DC') & !is.na(df$compoType) & df$owner == 'private', ] # 30701
  prop4 <- round(nrow(df4) * 0.04) # 1228
  df4 <- df4 %>% filter(structure == 'even', protect == 0, access == 1) %>% arrange(-rdi) # 3085 (8315 avec gini = 0.5)

  # what about 'C with with fir and/or spruce' and 'C' ?


  ###############################################################
  # save
  ###############################################################

  # reorder and rename columns
  df <- df[, c('cellID100', 'owner', 'access', 'protect', 'forestCellsPerHa',
               'compoType', 'gini', 'rdi', 'BA', 'Dg', 'meanH',	'structure',
               'density', 'stand')]
  #
  write.csv(df, paste0(landPath, '/managTable.csv'), row.names = F)

}
