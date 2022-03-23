managTable <- function(landscape){

  ###############################################################
  # initialisation
  ###############################################################

  # load packages
  require(reldist)
  require(dplyr)
  require(ggplot2)
  require(stringr)
  require(terra)

  # load environmental data
  env <- read.csv(paste0(landPath, '/cell25.csv'))

  # load virtual tree data
  tree <- read.csv(paste0(landPath, '/trees.csv'))
  # add cellID100 variable to tree data
  tree <- merge(tree, env[, c('cellID25', 'cellID100')], by = 'cellID25', all.x = TRUE, all.y = FALSE)
  tree <- tree[, c('cellID25', 'cellID100', 'sp', 'n', 'dbh', 'h')]

  if(landscape == 'bauges'){
    # load protected areas
    rb <- vect('./data/bauges/GEO/reserves_biologiques.shp')
    rn <- vect('./data/bauges/GEO/reserves_naturelles.shp')
    # load park limits
    park <- vect('./data/bauges/GEO/park.shp')
    # load ownership
    own <- vect('./data/bauges/GEO/Foret_publique_dep73-74_2814_dissolve.shp')
    # load accessibility
    access <- rast('./data/bauges/GEO/PNRfilled_F.distance.tif')
  } else if(landscape == 'milicz'){
    # load protected areas
    protect <- vect('./data/milicz/GEO/Milicz_protected_area.shp')
  } else if(landscape == 'sneznik'){
    elev <- rast(paste0(landPath, '/elev.asc'))
    managType <- vect('./data/sneznik/GEO/Sneznik_forest_stands2017_species_type.shp')
  }

  # load cellID100 raster
  cellID100 <- rast(paste0(landPath, '/cellID100.asc'))

  # load species self-thinning boundary parameters
  rdiParam <- read.table('./data/valeursCoefficientsRdi.txt', header = T)

  # load list of deciduous and coniferous sp
  deciduousSp <- readRDS('./data/deciduousSp.rds')
  coniferousSp <- readRDS('./data/coniferousSp.rds')


  ###############################################################
  # define protected areas
  ###############################################################

  if(landscape == 'bauges'){
    # keep only 'réserve biologique intégrale'
    # and exclude 'réserve biologique dirigée'
    rb <- rb[rb$code_r_enp == 'I',]
    # merge protected areas
    protect <- terra::union(rb, rn)
    protect$protection <- 1
    # select protected areas only in study area
    protect <- terra::intersect(protect, park) # TODO: different crs...

  }

  if(landscape == 'bauges' | landscape == 'milicz' ){
    # convert into raster
    # use cover to get proportion of each 100*100m cell covered by polygon
    protect <- rasterize(protect, cellID100, cover = TRUE, background = 0)
    names(protect) <- 'protect'

    # stack with cellID100
    protect <- c(cellID100, protect)

    # convert into dataframe
    protect <- as.data.frame(protect)

    # if protect >= 0.5 then most of the cell is covered by protected
    # area --> replace by 1. if protect < 0.5 --> replace by 0.
    df <- protect %>% mutate(protect = if_else(protect >= 0.5, 1, 0))
  }

  if(landscape == 'sneznik'){
    df <- as.data.frame(cellID100) %>% mutate(protect = 0)
  }

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

  if(landscape == 'bauges' | landscape == 'sneznik'){
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
  }

  if(landscape == 'bauges'){
  # identify deciduous and coniferous stands
  mainSp <- mainSp %>% group_by(cellID100) %>% mutate(D = sum(str_detect(mainSp, deciduousSp)),
                                                      C = sum(str_detect(mainSp, coniferousSp)),
                                                      compoType = if_else(D>0 & C>0, 'DC', if_else(D>0 & C==0, 'D', 'C')))
  #
  # classify main species into composition types:
    # beech
    mainSp[mainSp$mainSp == 'Fagus sylvatica' , 'compoType'] <- 'beech'
    # fir and or spruce
    fs <- c('Abies alba', 'Picea abies', 'Abies alba - Picea abies', 'Picea abies - Abies alba')
    mainSp[mainSp$mainSp %in% fs , 'compoType'] <- 'fir and or spruce'
    # mixed beech - fir and or spruce
    m <- c('Fagus sylvatica - Abies alba','Abies alba - Fagus sylvatica',
           'Fagus sylvatica - Picea abies','Picea abies - Fagus sylvatica',
           'Fagus sylvatica - Abies alba - Picea abies',
           'Fagus sylvatica - Picea abies - Abies alba',
           'Abies alba - Picea abies - Fagus sylvatica',
           'Abies alba - Fagus sylvatica - Picea abies',
           'Picea abies - Fagus sylvatica - Abies alba',
           'Picea abies - Abies alba - Fagus sylvatica')
    mainSp[mainSp$mainSp %in% m , 'compoType'] <- 'beech with fir and or spruce'

    # subdivide DC into DC and DC with fir and or spruce
    mainSp <- mainSp %>% mutate(DCfs = sum(str_detect(mainSp, c('Picea abies', 'Abies alba'))),
                                DCfsc = if_else(C > DCfs, 1, 0))
    mainSp[mainSp$compoType == 'DC' & mainSp$DCfs > 0 & mainSp$DCfsc == 0, 'compoType'] <- 'D with fir and or spruce'
    mainSp[mainSp$compoType == 'DC' & mainSp$DCfs > 0 & mainSp$DCfsc == 1, 'compoType'] <- 'DC with fir and or spruce'

    # subdivide C into C and C with fir and or spruce
    mainSp[mainSp$compoType == 'C' & mainSp$DCfs > 0, 'compoType'] <- 'C with fir and or spruce'

    # --------- group compoType based on the type of management they will follow
    # DC / C with fir and or spruce / DC with fir and or spruce / D with fir and or spruce ------> fir and or spruce with DC
    mainSp[mainSp$compoType %in% c('DC', 'C with fir and or spruce', 'DC with fir and or spruce', 'D with fir and or spruce'), 'compoType'] <- 'fir and or spruce with DC'
    # beech ------> D
    mainSp[mainSp$compoType == 'beech', 'compoType'] <- 'D'
    # --------- take into account other compo
    mainSp[!(mainSp$compoType %in% c('beech with fir and or spruce', 'D', 'fir and or spruce', 'fir and or spruce with DC')), 'compoType'] <- 'other compo'

  }

  if(landscape == 'milicz'){
    # retrieve dominant species
    mainSp <- tree %>% group_by(cellID100, sp) %>%
                              summarise(BA = sum((pi * (dbh/200)^2) * n)) %>%
                              group_by(cellID100) %>% arrange(cellID100, -BA) %>%
                              filter(BA == max(BA)) %>% rename(compoType = sp)
    # list of managed species
    spList <- c('Pinus sylvestris', 'Picea abies', 'Larix decidua',
                'Quercus robur', 'Fagus sylvatica', 'Betula pendula',
                'Alnus glutinosa')
    # if not managed species --> 'other sp'
    mainSp <- mainSp %>% mutate(compoType = case_when(compoType %in% spList ~ compoType, !(compoType %in% spList) ~ 'other sp'))

  }

  if(landscape == 'sneznik'){
    # add elevation
    elev100 <- aggregate(elev, 4, mean)
    elev100 <- c(cellID100, elev100)
    elev100 <- as.data.frame(elev100)
    mainSp <- left_join(mainSp, elev100, by = 'cellID100')
    mainSp <- mainSp %>% mutate(elev = case_when(elev >= 1100 ~ 'HA', elev < 1100 ~ 'LA'))

    # create composition types
    mainSp$compoType <- NA

    # beech dominated stands (beech is the most abundant species) + elevation
    mainSp[str_starts(mainSp$mainSp, 'Fagus sylvatica'), 'compoType'] <- 'beech dominated'
    mainSp[str_starts(mainSp$mainSp, 'Fagus sylvatica'), 'compoType'] <- apply(mainSp[str_starts(mainSp$mainSp, 'Fagus sylvatica'), c('compoType', 'elev')] , 1 , paste , collapse = "-" )

    # fir or spruce dominated stands (fir OR spruce most abundant species)
    mainSp[str_starts(mainSp$mainSp, 'Abies alba'), 'compoType'] <- 'fir or spruce dominated'
    mainSp[str_starts(mainSp$mainSp, 'Picea abies'), 'compoType'] <- 'fir or spruce dominated'

    # fir-beech-spruce mixed stands
    mainSp[str_detect(mainSp$mainSp, 'Abies alba'), 'fir'] <- 1
    mainSp[str_detect(mainSp$mainSp, 'Fagus sylvatica'), 'beech'] <- 1
    mainSp[str_detect(mainSp$mainSp, 'Picea abies'), 'spruce'] <- 1
    # # fir AND beech are most abundant species
    # mainSp$mixfb <- apply(mainSp[, c('fir', 'beech')], 1, sum, na.rm = TRUE)
    # mainSp[mainSp$mixfb == 2, 'compoType'] <- 'mixed'
    # # spruce AND beech are most abundant species
    # mainSp$mixsb <- apply(mainSp[, c('beech', 'spruce')], 1, sum, na.rm = TRUE)
    # mainSp[mainSp$mixsb == 2, 'compoType'] <- 'mixed'
    # fir AND beech AND spruce are most abundant species
    mainSp$mixfbs <- apply(mainSp[, c('fir', 'beech', 'spruce')], 1, sum, na.rm = TRUE)
    mainSp[mainSp$mixfbs == 3, 'compoType'] <- 'mixed'

  }

  # add to df
  df <- merge(df, mainSp[, c('cellID100', 'compoType')], by = 'cellID100', all.x = TRUE)


  ###############################################################
  # ownership
  ###############################################################

  # retrieve ownership in study area
  if(landscape == 'bauges'){
    # convert into raster
    # use cover to get proportion of each 100*100m cell covered by polygon
    own <- crop(own, park)
    own <- terra::rasterize(own, cellID100, cover = TRUE, background = 0)
    names(own) <- 'public'

    # stack with cellID100
    own <- c(cellID100, own)

    # convert into dataframe
    own <- terra::as.data.frame(own)

    # if public >= 0.5 then most of the cell is public --> replace by public.
    # if public < 0.5 --> replace by private.
    own <- own %>% mutate(owner = if_else(public >= 0.5, 'public', 'private'))

    # add to df
    df <- left_join(df, own[, c('cellID100', 'owner')], by = 'cellID100')

  }

  if(landscape == 'milicz'){
    # we only consider public forests in milicz as defined in saveLandscape
    # --> set ownership to public wherever there's forest
    df <- df %>% mutate(owner = ifelse(!is.na(BA), 'public', NA))

  }

  if(landscape == 'sneznik'){
    # management does not depend on ownership
    # --> set ownership to NA
    df$owner <- NA
  }


  ###############################################################
  # access
  ###############################################################

  if(landscape == 'bauges'){

    # set shp extent
    access <- crop(access, cellID100)

    # aggregate access value from 5*5m to 100*100m
    access <- terra::aggregate(access, fact = 20, fun = 'mean', na.rm = TRUE)

    # stack with cellID100
    access <- c(cellID100, access)

    # convert into dataframe
    access <- terra::as.data.frame(access)
    names(access) <- c('cellID100', 'dist')

    # if access > 2km --> not accessible
    access <- access %>% mutate(access = if_else(dist <= 2000, 1, 0)) %>%
                         dplyr::select(-dist)
    # access = 1 --> accessible
    # access = 0 --> inaccessible

    # add to df
    df <- left_join(df, access, by = 'cellID100')
    df[is.na(df$access), 'access'] <- 0

  }

  if(landscape == 'milicz' | landscape == 'sneznik'){
    df$access <- 1
  }

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

  ###############################################################
  # define density class
  ###############################################################

  # calculate BA/ha
  df <- df %>% mutate(BA_ha = 16 * BA / forestCellsPerHa)

  if(landscape == 'bauges'){
    # uneven-aged stands
    # threshold of BA after thinning from the sylviculture guide for mountain forest
    # conifers: 35 - 30 - 25 m2 for high - medium - low density
    # mixed: 30 - 25 - 20 m2
    # deciduous: 25 - 20 - 15 m2

    # conifers
    df[df$structure == 'uneven' & !is.na(df$structure) & !is.na(df$compoType) & df$compoType == 'fir and or spruce' & df$BA_ha < 30, 'density'] <- 'low'
    df[df$structure == 'uneven' & !is.na(df$structure) & !is.na(df$compoType) & df$compoType == 'fir and or spruce' & df$BA_ha >= 30 & df$BA_ha < 35, 'density'] <- 'medium'
    df[df$structure == 'uneven' & !is.na(df$structure) & !is.na(df$compoType) & df$compoType == 'fir and or spruce' & df$BA_ha >= 35, 'density'] <- 'high'

    # mixed
    comp <- c('fir and or spruce with DC', 'beech with fir and or spruce')
    df[df$structure == 'uneven' & !is.na(df$structure) & !is.na(df$compoType) & df$compoType %in% comp & df$BA_ha < 25, 'density'] <- 'low'
    df[df$structure == 'uneven' & !is.na(df$structure) & !is.na(df$compoType) & df$compoType %in% comp & df$BA_ha >= 25 & df$BA_ha < 30, 'density'] <- 'medium'
    df[df$structure == 'uneven' & !is.na(df$structure) & !is.na(df$compoType) & df$compoType %in% comp & df$BA_ha >= 30, 'density'] <- 'high'

    # deciduous
    df[df$structure == 'uneven' & !is.na(df$structure) & !is.na(df$compoType) & df$compoType == 'D' & df$BA_ha < 20, 'density'] <- 'low'
    df[df$structure == 'uneven' & !is.na(df$structure) & !is.na(df$compoType) & df$compoType == 'D' & df$BA_ha >= 20 & df$BA_ha < 25, 'density'] <- 'medium'
    df[df$structure == 'uneven' & !is.na(df$structure) & !is.na(df$compoType) & df$compoType == 'D' & df$BA_ha >= 25, 'density'] <- 'high'

    # even-aged stands
    # mean target rdi: 0.6 and 0.7
    # after thinning targets: 0.55 and 0.65
    # stands with rdi up to 0.65 are considered to have a target rdi of 0.6
    # stands with rdi > 0.65 are considered to have a target rdi of 0.7
    df[df$structure == 'even' & !is.na(df$structure) & df$rdi < 0.65, 'density'] <- 'medium'
    df[df$structure == 'even' & !is.na(df$structure) & df$rdi >= 0.65, 'density'] <- 'high'

    #---------------------------------------------------
    # add 'final cut' management based on BA, rdi, Dg
    df$manag <- paste(df$structure, '-', df$density)

    # uneven-aged stands
    # coniferous stands above 50 m2 are considered abandoned and will only undergo a final cut
    df[df$structure == 'uneven' & !is.na(df$structure) & !is.na(df$compoType) & df$compoType == 'fir and or spruce' & df$BA_ha >= 50, 'manag'] <- 'final cut'
    # mixed stands above 45 m2 are considered abandoned and will only undergo a final cut
    df[df$structure == 'uneven' & !is.na(df$structure) & !is.na(df$compoType) & df$compoType %in% comp & df$BA_ha >= 45, 'manag'] <- 'final cut'
    # mixed stands above 45 m2 are considered abandoned and will only undergo a final cut
    df[df$structure == 'uneven' & !is.na(df$structure) & !is.na(df$compoType) & df$compoType == 'D' & df$BA_ha >= 40, 'manag'] <- 'final cut'

    #  even-aged stands
    # dense stands (rdi > 0.8) with big trees (Dg > 30) are considered abandoned and will only undergo a final cut
    df[df$structure == 'even' & !is.na(df$structure) & df$rdi >= 0.8 & df$Dg >= 20, 'manag'] <- 'final cut'
    # XX% of deciduous stands with small trees (Dg < 20) are concidered coppice (whatever their rdi)
    # create column with integer for subsetting those stands
    df$sub <- rep(x = c(1,2), length.out = nrow(df)) # change x to change size of subset
    df[df$structure == 'even' & !is.na(df$structure) & df$Dg <= 20 & df$compoType == 'D' & !is.na(df$compoType) & df$sub == 2, 'manag'] <- 'coppice'

    #---------------------------------------------------
    # other compo -> no manag & inaccessible
    df[!is.na(df$compoType) & df$compoType == 'other compo', 'manag'] <- 'no manag'
    df[!is.na(df$compoType) & df$compoType == 'other compo', 'access'] <- 0

  }

  if(landscape == 'milicz'){
    df <- df %>% mutate(density = NA,
                        manag = case_when(compoType != 'other sp'~ compoType, compoType == 'other sp'~ 'no manag'))

  }

  if(landscape == 'sneznik'){
    df <- df %>% mutate(density = NA)
    # assign management type
    manag <- raster::rasterize(managType, cellID100, fun = function(x,...) max(x, na.rm = TRUE), field = 'StandType')
    names(manag) <- 'manag'
    manag <- c(cellID100, manag)
    manag <- as.data.frame(manag)
    # if manag > 1 --> uneven
    # if manag <=1 --> even
    manag <- manag %>% mutate(manag = case_when(manag > 1 ~ 'uneven', manag <= 1 ~ 'even'))
    df <- left_join(df, manag, by = 'cellID100')
    # assign most common manag type to cells with trees but manag = NA
    df[!is.na(df$BA) & is.na(df$manag), 'manag'] <- 'even'
    # remove manag type for cells with no trees
    df[is.na(df$BA), 'manag'] <- NA
    # assign only manag type 'even' for compoType == 'beech dominated'
    types <- c('beech dominated-LA', 'beech dominated-HA')
    df[df$compoType %in% types & !is.na(df$compoType), 'manag'] <- 'even'

  }

  ###############################################################
  # define 'no management' areas
  ###############################################################

  df[(df$access == 0 | df$protect == 1) & !is.na(df$forestCellsPerHa) & !is.na(df$access) & !is.na(df$protect), 'manag'] <- 'no manag'

  ###############################################################
  # save
  ###############################################################

  # reorder columns
  colOrd <- c('cellID100', 'owner', 'access', 'protect', 'forestCellsPerHa',
               'compoType', 'gini', 'rdi', 'BA', 'BA_ha', 'Dg', 'meanH',	'structure',
               'density', 'manag')
  # round numbers
  colRou <- c('gini', 'rdi', 'BA', 'BA_ha', 'Dg', 'meanH')
  df <- df %>% dplyr::select(all_of(colOrd)) %>% arrange(cellID100) %>%
               mutate(across(all_of(colRou), round, 4))
  #
  write.csv(df, paste0(landPath, '/managTableCell100.csv'), row.names = F)

}
