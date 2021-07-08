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
library(tidyr)

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

###############################################################
# define whether the slope makes logging impossible
###############################################################

# calculate mean slope on 1ha sites
# logging impossible if slope >= 110% (slope > 47.73 degrees)
df <- env %>% group_by(cellID100) %>% summarise(slope = mean(slope)) %>%
              mutate(loggable = if_else(slope < 47.73, 1, 0)) %>%
              dplyr::select(-slope)
#
# loggable = 1 --> can be logged
# loggable = 0 --> cannot be logged

###############################################################
# calculate Gini index on 1ha cells and define stand type
# (evenaged / unevenaged)
###############################################################

# if gini<0.45 --> evenaged stand, else --> unevenaged stand
gini <- tree %>% group_by(cellID100) %>% mutate(ba = pi * dbh^2 / 4) %>%
                 summarise(gini = gini(x = ba, weights = n),
                           BA_m2 = sum((pi * (dbh/200)^2) * n),
                           Dg_cm = sqrt(sum(dbh^2 * n)/sum(n)),
                           meanH_m = sum(h * n) /sum(n)) %>%
                 mutate(type = if_else(gini < 0.45, 'even', 'uneven')) %>%
                 dplyr::select(-gini)
#

# add to df
df <- merge(df, gini, by = 'cellID100', all.x = TRUE)
df$type <- as.factor(df$type)

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
protect <- union(rb, rn)
protect$protection <- 1

# select protected areas only in study area
protect <- intersect(protect, park)

# convert into raster
# use getCover to define proportion of each 100*100m cell covered by polygon
protect <- rasterize(protect, cellID100, getCover = TRUE)
names(protect) <- 'protect'

# stack with cellID100
protect <- stack(cellID100, protect)

# convert into dataframe
protect <- as.data.frame(protect)

# if protect >= 0.5 then most of the cell is covered by protected
# area --> replace by 1. if protect < 0.5 --> replace by 0.
protect <- protect %>% mutate(protect = if_else(protect >= 0.5, 1, 0))

# add to df
df <- merge(df, protect, by = 'cellID100')


###############################################################
# define composition type
###############################################################

# set threshold to identify mainSp
# thresh = 0.8 means you will get the species making up for >= 80% of
# the stand BA
thresh = 0.70
# retrieve main species
mainSp <- tree %>% group_by(cellID100, sp) %>%
                          summarise(BA = sum((pi * (dbh/200)^2) * n)) %>%
                          group_by(cellID100) %>% arrange(cellID100, -BA) %>%
                          mutate(BAtot = sum(BA), BAprop = BA/BAtot, cumulProp = cumsum(BAprop), HigherThanthresh = case_when(
                            cumulProp >= thresh ~ cumulProp), minCompo = min(HigherThanthresh, na.rm = TRUE)) %>%
                          filter(cumulProp <= minCompo) %>% summarise(mainSp = paste(sp, collapse=' - '))
#

# class main species into composition types:
# deciduous and coniferous sp list----------------------------------------------
deciduousSp <- c('Fagus sylvatica', 'Sorbus aria', 'Quercus robur',
                 'Quercus petraea', 'Acer opalus', 'Tilia platyphyllos',
                 'Carpinus betulus', 'Castanea sativa', 'Populus nigra',
                 'Fraxinus excelsior', 'Corylus avellana', 'Acer campestre',
                 'Betula pendula', 'Populus tremula', 'Prunus avium',
                 'Acer pseudoplatanus', 'Tilia cordata', 'Salix caprea',
                 'Ulmus glabra', 'Sorbus aucuparia', 'Robinia pseudacacia',
                 'Malus sylvestris', 'Acer platanoides', 'Salix alba',
                 'Crataegus monogyna', 'Sorbus mougeotii', 'Laburnum anagyroides',
                 'Quercus pubescens', 'Alnus glutinosa', 'Salix cinerea',
                 'Alnus incana', 'Ilex aquifolium', 'Buxus sempervirens',
                 'Prunus padus', 'Robinia pseudoacacia')
#
coniferousSp <- c('Picea abies', 'Abies alba', 'Pinus sylvestris',
                  'Taxus baccata', 'Larix decidua', 'Larix kaempferi',
                  'Picea sitchensis')
#-------------------------------------------------------------------------------
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

# add to df
df <- merge(df, mainSp[, c('cellID100', 'compoType')], by = 'cellID100', all.x = TRUE)

# plot composition types proportion in landscape
compo <- mainSp %>% group_by(compoType) %>% mutate(surface = 1) %>%
                    summarise(surface = sum(surface)) %>% arrange(-surface) %>%
                    ungroup() %>% mutate(totSurf = sum(surface)) %>%
                    group_by(compoType) %>% mutate(relSurf = surface * 100 / totSurf) %>%
                    dplyr::select(-totSurf)
compo$compoType <- factor(compo$compoType, levels = compo$compoType)

# plot composition at 1ha cells
ggplot(data = compo) +
geom_bar(aes(x = compoType, y = surface), stat = 'identity') +
geom_text(aes(x = compoType, y = surface, label = paste(round(relSurf, 2), '%') ), vjust = -0.3, size = 5, col = 'orange') +
geom_text(aes(x = compoType, y = surface, label = paste(round(surface, 2), 'ha') ), vjust = +1.2, size = 5, col = 'orange') +
theme_light() +
xlab('main species') +
ylab('surface (ha)') +
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
     plot.title = element_text(hjust = 0.5))
#

# TODO: sortir liste des espèces decidues/coniferes (ici et dans spTransform)

###############################################################
# ownership
###############################################################

# retrieve ownership in study area
own <- intersect(own, park)

# convert into raster
# use getCover to define proportion of each 100*100m cell covered by polygon
own <- rasterize(own, cellID100, getCover = TRUE)
names(own) <- 'public'

# stack with cellID100
own <- stack(cellID100, own)

# convert into dataframe
own <- as.data.frame(own)

# if public >= 0.5 then most of the cell is public --> replace by 1.
# if public < 0.5 --> replace by 0.
own <- own %>% mutate(public = if_else(public >= 0.5, 1, 0))

# public = 1 --> public
# public = 0 --> private

# add to df
df <- merge(df, own, by = 'cellID100')


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

# TODO: clarifier overlap entre loggable et access
# table(df$loggable, df$access)
# hist(df[df$loggable == 0, 'access'])
# hist(df[df$loggable == 1, 'access'])


###############################################################
# number of forest 25*25m cells in 100*100m cells
###############################################################

# nb of forest cells per ha
forCel <- tree %>% group_by(cellID100) %>% summarise(forCel = length(unique(cellID25)))

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
rdi <- merge(rdi, df[, c('cellID100', 'forCel')], by = 'cellID100', all.x = TRUE)

# calculate species Nmax
rdi$Nspmax <- exp( rdi$rqIntercept + rdi$rqSlope * log(rdi$Dgsp) ) * (rdi$forCel / 16)

# Calculate sp partial rdi
rdi$rdip <- rdi$nsp / rdi$Nspmax

# claculate total rdi
rdit <- rdi %>% group_by(cellID100) %>% summarise(rdi = sum(rdip))

# add to df
df <- merge(df, rdit, by = 'cellID100', all.x = TRUE)

# define density class for uneven-aged stands
qtuneven <- quantile(df[df$type == 'uneven' & df$access == 1, 'rdi'], na.rm = T, c(0.33, 0.66))
hist(df[df$type == 'uneven' & df$access == 1, 'rdi'], breaks = 100)
abline(v = qtuneven, col = 'red', lty = 5, lwd = 2)
df[df$type == 'uneven' & !is.na(df$type), 'density'] <- 'medium'
df[df$type == 'uneven' & !is.na(df$type) & df$rdi > qtuneven[2], 'density'] <- 'high'
df[df$type == 'uneven' & !is.na(df$type) & df$rdi < qtuneven[1], 'density'] <- 'low'

# define density class for even-aged stands
qteven <- quantile(df[df$type == 'even' & df$access == 1, 'rdi'], na.rm = T, 0.5)
hist(df[df$type == 'even' & df$access == 1, 'rdi'], breaks = 100)
abline(v = qteven, col = 'red', lty = 5, lwd = 2)
df[df$type == 'even' & !is.na(df$type), 'density'] <- 'low'
df[df$type == 'even' & !is.na(df$type) & df$rdi > qteven, 'density'] <- 'high'

# TODO: corriger rdi là où il n'y a pas 16 cellules de foret

###############################################################
# define stand type (compoType + even / uneven / protect / ...)
###############################################################

# no management
noMan <- df %>% filter(!is.na(compoType)) %>% dplyr::select(cellID100, compoType, access, public) %>%
                filter(access == 0) %>%
                pivot_wider(names_from = compoType, values_from = access) %>%
                dplyr::select(-cellID100)
noManPub <- noMan %>% filter(public == 1) %>% dplyr::select(-public) %>%
                      summarise_all(~sum(!is.na(.)))
noManPri <- noMan %>% filter(public == 0) %>% dplyr::select(-public) %>%
                      summarise_all(~sum(!is.na(.)))
#



--> ajouter distinction public / privé / densité

# uneven-aged
uneven <- df %>% filter(!is.na(compoType)) %>% dplyr::select(cellID100, compoType, type) %>%
                filter(type == 'uneven') %>%
                pivot_wider(names_from = compoType, values_from = type) %>%
                dplyr::select(-cellID100) %>%
                summarise_all(~sum(!is.na(.)))
#
# even-aged
even <- df %>% filter(!is.na(compoType)) %>% dplyr::select(cellID100, compoType, type) %>%
                filter(type == 'even') %>%
                pivot_wider(names_from = compoType, values_from = type) %>%
                dplyr::select(-cellID100) %>%
                summarise_all(~sum(!is.na(.)))
#
# public

#
ggplot(data = df[!is.na(df$compoType),]) +
geom_bar(aes(x = compoType), stat = 'count')



# # assign compoType to stand type
# df$standType <- paste(df$compoType, df$type)
#
# # mark as 'nologging' non-loggable stands
# df[df$loggable == 0 & !is.na(df$loggable), 'standType'] <- paste(df[df$loggable == 0 & !is.na(df$loggable), 'compoType'], 'noLogging')
#
# # mark as 'protected' protected stands
# df[df$protect == 1, 'standType'] <- 'protected'
#
# # plot stand type share in the landscape
# typeShare <- df %>% mutate(surface = 1) %>% group_by(standType) %>% filter(compoType != 'NA') %>%
#                     summarise(surface = sum(surface)) %>% arrange(-surface) %>%
#                     ungroup() %>% mutate(totSurf = sum(surface)) %>%
#                     group_by(standType) %>% mutate(relSurf = surface * 100 / totSurf) %>%
#                     dplyr::select(-totSurf)
#
# typeShare$standType <- factor(typeShare$standType, levels = typeShare$standType)
# #
# ggplot(data = typeShare) +
# geom_bar(aes(x = standType, y = surface), stat = 'identity') +
# geom_text(aes(x = standType, y = surface, label = paste(round(relSurf, 2), '%') ), vjust = -0.3, size = 5, col = 'orange') +
# geom_text(aes(x = standType, y = surface, label = paste(round(surface, 2), 'ha') ), vjust = +1.2, size = 5, col = 'orange') +
# theme_light() +
# xlab('main species') +
# ylab('surface (ha)') +
# theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
#      plot.title = element_text(hjust = 0.5))
# #
#
# # plot on map
# cellID100$standType <- as.numeric(as.factor(df$standType))
# plot(cellID100$standType)
# plot(park, add = T)




# TODO: vérifier que toutes les forets ont des valeurs de access/protect...
# croiser les fino des colonnes dans les deux sens e.g. access->foret, foret->access
