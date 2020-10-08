###############################################################
# initialisation
###############################################################

# clean up environment
rm(list = ls())

# load packages
library(rgdal)
library(raster)
library(dplyr)
library(ggplot2)

# set work directory
setwd("C:/Users/raphael.aussenac/Documents/GitHub/LandscapeInit")

# load forest cover data
forest <- raster("./Init/forest.asc")

# load TFV spatial data
bd <- readOGR(dsn = "./data/GEO", layer = "BD_Foret_V2_PNRfilled_Foret_2014", encoding = "UTF-8", use_iconv = TRUE)

# load NFI tree data
tree <- read.csv('./data/NFI/BaugesTrees.csv', sep = ';')

# load correspondence between NFI plots and TFV types
crpdTFV <- read.csv('./data/NFI/croisementIfn.csv', sep = ';')

###############################################################
# group TFV types together
###############################################################

# first assign TFV to each NFI plot using crpdTFV
tree <- merge(tree, crpdTFV[, c('idp', 'TFV')], by.x = 'id_plot', by.y = 'idp')
tree$TFV <- droplevels(tree$TFV)
colnames(tree)[colnames(tree) == 'TFV'] <- 'CODE_TFV'

# Surface area of each TFV type and number of associated NFI plots
TFVcountAndSurface <- function(bd, tree){
  # calculate area of each TFV type
  bd$area <- area(bd)
  surfArea <- data.frame(bd) %>% group_by(CODE_TFV) %>%
                                summarize(surface = sum(area)/10000) %>%
                                arrange(-surface)
  # count number of NFI in each TFV type
  TFVplotCount <- tree %>% group_by(CODE_TFV) %>% summarise(N = length(unique(id_plot)))
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
  geom_text(data = surfArea[3,], aes(x = CODE_TFV, y = surface, label =  paste(round(prop, 2), "% of forest cover")), vjust = -1, hjust = 0, size = 3.5, col = 'black') +
  # annotate("text", x = 3, y = 10000, label = paste(round(prop, 2), "% of forest cover")) +
  theme_light() +
  xlab('TFV code') +
  ylab('surface (ha)') +
  # scale_y_continuous("surface (ha)", sec.axis = sec_axis(~ . * 1.20, name = "number of NFI plots")) +
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
# not specified --------------------------------------------> A

groupTfv <- function(df){
  df <- df[!(df$CODE_TFV %in% c('FO0', 'FF0', 'LA4', 'LA6')),]
  df[df$CODE_TFV == "FO1", "CODE_TFV"] <- A
  df[df$CODE_TFV == "FF1-00", "CODE_TFV"] <- A
  df[df$CODE_TFV == "FO2", "CODE_TFV"] <- B
  df[df$CODE_TFV == "FO3", "CODE_TFV"] <- C
  df[df$CODE_TFV == "FF1-49-49", "CODE_TFV"] <- A
  df[df$CODE_TFV == "FF2-00-00", "CODE_TFV"] <- B
  df[df$CODE_TFV == "FF1G01-01", "CODE_TFV"] <- A
  df[df$CODE_TFV == "FF2G53-53", "CODE_TFV"] <- B
  df[df$CODE_TFV == "FF2-90-90", "CODE_TFV"] <- B
  df[df$CODE_TFV == "FF2-64-64", "CODE_TFV"] <- B
  df[df$CODE_TFV == "FF2-00", "CODE_TFV"] <- B
  df[df$CODE_TFV == "FF1-10-10", "CODE_TFV"] <- A
  df[df$CODE_TFV == "FF2-52-52", "CODE_TFV"] <- B
  df[df$CODE_TFV == "", "CODE_TFV"] <- A
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











# convert forest into raster and set extent + resolution
ext <- floor(extent(forest))
r <- raster(ext, res=res(forest))
bd$CODE_TFV <- as.factor(bd$CODE_TFV)
tfvRaster <- rasterize(bd, r, field = "CODE_TFV")
# set projection
crs(tfvRaster) <- crs(forest)


###############################################################
# assign new simplified TFV to all NFI plots
###############################################################



###############################################################
# assign a composition to each cell
###############################################################

# define coniferous and/or deciduous dominant species at each site
