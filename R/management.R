###############################################################
# initialisation
###############################################################

# load packages
require(reldist)
require(dplyr)
require(rgdal)
require(raster)

# load virtual tree data
tree <- read.csv(paste0(landPath, '/trees75.csv'))

# load environmental data
env <- read.csv(paste0(landPath, '/envVariables.csv'))

# load protected areas
rb <- readOGR(dsn = './data/bauges/GEO', layer = 'reserves_biologiques', encoding = 'UTF-8', use_iconv = TRUE)
rn <- readOGR(dsn = './data/bauges/GEO', layer = 'reserves_naturelles', encoding = 'UTF-8', use_iconv = TRUE)
park <- readOGR(dsn = './data/bauges/GEO', layer = 'park', encoding = 'UTF-8', use_iconv = TRUE)

# load cellID100 raster
cellID100 <- raster(paste0(landPath, '/cellID100.asc'))

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
# use getCover to define proportion of each cell covered by polygon
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

--> prendre script composition dans evalcompo
- beech
- fir and/or spruce
- mixed beech - fir and/or spruce
- other deciduous
---- ajouter oak?

###############################################################
# number of forest 25*25m cells in ha
###############################################################



###############################################################
# calculate rdi
###############################################################

valeurs sp / feuillus / pins / autres résineux



###############################################################
###############################################################
###############################################################

# lab

hist(protect[protect$protect > 0,])



plot(rb, col = 'green', border = 'green')
plot(ri, col = 'red', border = 'red', add =TRUE)
plot(rn, col = 'blue', border = 'blue', add =TRUE)
plot(cellID25, add = TRUE)


irr <- gini %>% filter(Gini > 0.45)
reg <- gini %>% filter(Gini <= 0.45)


nrow(irr) * 100 /nrow(gini)

83% d'irrégulier
16% de regulier



hist(df$Gini, breaks = 100)
range(df$Gini)



# library(plotrix)
weighted.hist(x = test$dbh, w = test$n)


8636 -> 0.270


test <- tree[tree$cellID100 == 8636, ]
gini(x = test$dbh, weights = test$n)
