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
library(plotly)
library(doParallel)

# set work directory
setwd("C:/Users/raphael.aussenac/Documents/GitHub/LandscapeInit")

# load LIDAR rasters
Dg <- raster("./Init/Dg.asc")
BA <- raster("./Init/BA.asc")
Dprop <- raster("./Init/Dprop.asc")
Dprop <- Dprop / 100

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


###############################################################
# Retrieve Dg, BA, Dprop and TFV for all forest cells
###############################################################

# first assign TFV code to each forest cell
bd$CODE_TFV <- factor(bd$CODE_TFV, levels = c('FF1-00-00', 'FF2G61-61', 'FF31', 'FF32', 'FF1-09-09'))
bd$CODE_TFV <- as.numeric(bd$CODE_TFV)
TFVraster <- rasterize(bd, Dg, field = "CODE_TFV")
# merge with deciduous proportion of BA
compoRaster <- stack(TFVraster, Dprop, Dg, BA)

plot(compoRaster$layer == 1 & compoRaster$Dprop >= 50 & compoRaster$dg > 20)

###############################################################
# Retrieve Dg, BA, Dprop and TFV for all NFI plots
###############################################################

# correct species name
tree$species_name <- as.character(tree$species_name)
tree[substr(tree$species_name, 1, 4) == "Ilex", 'species_name'] <- 'Ilex aquifolium'

# first define deciduous and coniferous species
deciduousSp <- c('Fagus sylvatica', 'Sorbus aria', 'Quercus robur',
                 'Quercus petraea', 'Acer opalus', 'Tilia platyphyllos',
                 'Carpinus betulus', 'Castanea sativa', 'Populus nigra',
                 'Fraxinus excelsior', 'Corylus avellana', 'Acer campestre',
                 'Betula pendula', 'Populus tremula', 'Prunus avium',
                 'Acer pseudoplatanus', 'Tilia cordata', 'Salix caprea',
                 'Ulmus glabra', 'Sorbus aucuparia', 'Robinia pseudacacia',
                 'Malus sylvestris', 'Acer platanoides', 'Salix alba',
                 'Crataegus monogyna', 'Sorbus mougeoti', 'Laburnum anagyroides',
                 'Quercus pubescens', 'Alnus glutinosa', 'Salix cinerea',
                 'Alnus incana', 'Ilex aquifolium', 'Buxus sempervirens',
                 'Prunus padus')
#
coniferousSp <- c('Picea abies', 'Abies alba', 'Pinus sylvestris',
                  'Taxus baccata', 'Larix decidua', 'Larix kaempferi')
#
tree[tree$species_name %in% deciduousSp, 'spType'] <- 'D'
tree[tree$species_name %in% coniferousSp, 'spType'] <- 'C'

# calculate proportion of deciduous BA
tree$DBH <- tree$c130_final / pi
NFIDprop <- tree %>% group_by(id_plot, spType) %>%
                          summarise(BAdc = sum((pi * (DBH/200)^2) * weight)) %>%
                          group_by(id_plot) %>% mutate(BA = sum(BAdc), Dprop = BAdc/BA) %>%
                          filter(spType == 'D') %>% select(-spType, -BAdc, -BA)
# calculate Dg and BA
NFI <- tree %>% group_by(id_plot) %>%
                    summarise(Dg = sqrt(sum(DBH^2 * weight)/sum(weight)),
                              BA = sum((pi * (DBH/200)^2) * weight),
                              CODE_TFV = unique(CODE_TFV))
NFI <- merge(NFIDprop, NFI, by = 'id_plot', all = TRUE)
NFI[is.na(NFI$Dprop), 'Dprop'] <- 0

###############################################################
# normalise Dg and BA values to 0-1 range
# range defined on NFI values
###############################################################

# NFI plots
NFI$Dg01 <- ( NFI$Dg - min(NFI$Dg) ) / ( max(NFI$Dg) - min(NFI$Dg) )
NFI$BA01 <- ( NFI$BA - min(NFI$BA) ) / ( max(NFI$BA) - min(NFI$BA) )

# plot 3d
fig <- plot_ly(x = NFI$Dg01, y = NFI$Dprop, z = NFI$BA01, type="scatter3d", mode = "markers", color = NFI$CODE_TFV)
axx <- list(title = "quadratic diameter (cm)")
axy <- list(title = "proportion of Deciduous basal area")
axz <- list(title = "basal area (m^2/ha)")
fig <- fig %>% layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
fig

# forest cells
compoRaster$Dg01 <- ( compoRaster$dg - min(NFI$Dg) ) / ( max(NFI$Dg) - min(NFI$Dg) )
compoRaster$BA01 <- ( compoRaster$BA - min(NFI$BA) ) / ( max(NFI$BA) - min(NFI$BA) )

###############################################################
# claculate distance between forest cells values and NFI values
# and assign composition of nearest NFI plot
###############################################################

# first create list of NFI plot for each TFV type
plotsInTFV <- tree %>% group_by(id_plot) %>% summarize(CODE_TFV = unique(CODE_TFV)) %>% arrange(CODE_TFV)
plotsInTFV$CODE_TFV <- factor(plotsInTFV$CODE_TFV, levels = c('FF1-00-00', 'FF2G61-61', 'FF31', 'FF32', 'FF1-09-09'))
plotsInTFV$CODE_TFV <- as.numeric(plotsInTFV$CODE_TFV)

# convert NFI TFC CODE in numeric values
NFI$CODE_TFV <- factor(NFI$CODE_TFV, levels = c('FF1-00-00', 'FF2G61-61', 'FF31', 'FF32', 'FF1-09-09'))
NFI$CODE_TFV <- as.numeric(NFI$CODE_TFV)
NFI <- NFI %>% select(id_plot, Dprop, Dg01, BA01, CODE_TFV)

# add extra layer to compoRaster to recieve compo values (= ref stand id)
compo <- BA
compo[!is.na(compo)] <- NA
names(compo) <- "compo"
compoRaster <- stack(compoRaster, compo)

################################################################################
# # simple for loop
# test <- crop(compoRaster, extent(compoRaster)/20)
#
# start <- Sys.time()
#
# for (i in 1:nrow(test[])){
#   # print(i)
#   # some cells do not have any TFV value
#   if(!is.na(test[i][1])){
#     # retrieve list of plots (and their Dg, BA, Dprop values) in the TFV
#     # type associated to the cell
#     plotlist <- NFI %>% filter(CODE_TFV == test[i][1]) %>% select(-CODE_TFV)
#     rownames(plotlist) <- plotlist$id_plot
#     plotlist <- plotlist %>% select(-id_plot)
#     # calculate distance
#     plotlist$DpropCell <- test[i][2]
#     plotlist$Dg01Cell <- test[i][5]
#     plotlist$BA01Cell <- test[i][6]
#     plotlist$distance <- apply(plotlist, 1, function(x) dist(matrix(x, nrow = 2, byrow = TRUE)))
#     NearestPlot <- as.numeric(rownames(plotlist[plotlist$distance == min(plotlist$distance),]))
#     # assign nearest plot value to cell
#     test[i][7] <- NearestPlot
#   }
# }
#
# end <- Sys.time()
# end - start
################################################################################
############################## parallel



miniBauges <- crop(compoRaster, extent(compoRaster)/1.5)

miniBauges <- compoRaster

library(SpaDES)
essai <- splitRaster(compoRaster, nx = 2, ny = 2, buffer = 0)#, path, cl)
SplitRas(raster=compoRaster, ppside=3,save=FALSE,plot=TRUE)

split_raster(
  "./Init/Dg.asc",
  s = 2,
  "./Init/",
  gdalinfoPath = NULL,
  gdal_translatePath = NULL
)


blabla <- function(cell, i){
  # print(i)
  # some cells do not have any TFV value
  if(!is.na(cell[1])){
    # retrieve list of plots (and their Dg, BA, Dprop values) in the TFV
    # type associated to the cell
    plotlist <- NFI[NFI$CODE_TFV == cell[1], ]
    rownames(plotlist) <- plotlist$id_plot
    plotlist[, c('id_plot', 'CODE_TFV')] <- NULL
    # calculate distance
    plotlist$DpropCell <- cell[2]
    plotlist$Dg01Cell <- cell[5]
    plotlist$BA01Cell <- cell[6]
    plotlist$distance <- apply(plotlist, 1, function(x) dist(matrix(x, nrow = 2, byrow = TRUE)))
    NearestPlot <- as.numeric(rownames(plotlist[plotlist$distance == min(plotlist$distance),]))
    NearestPlot <- NearestPlot[1]
  } else {
    NearestPlot <- NA
  }
  # assign nearest plot value (or NA) to cells
  return(c(i,NearestPlot))
}

# i <- 1
# for (i in 1:5){
#   print(blabla(cell = miniBauges[i], i = i))
# }

cl <- makeCluster(6)
registerDoParallel(cl)
n <- nrow(miniBauges[])
start <- Sys.time()
results <- foreach(i = 1:n, .combine = 'rbind', .packages = c('raster', 'rgdal')) %dopar% {blabla(cell = miniBauges[i], i = i)}
end <- Sys.time()
end - start
stopCluster(cl)


results <- data.frame(results)
colnames(results) <- c('i', 'id')
results <- results %>% arrange(i)
miniBauges$compo <- results[,2]
plot(miniBauges)


#               extent             10      5       4      2
temps <- data.frame('nCells' = c(32760, 132158, 206358, 826276, 1469430),
                    'temps' = c(15.96, 65.4, 100.6, 569.4, 1344.6))

plot(temps$nCells, temps$temps, type = 'l')
points(temps$nCells, temps$temps, pch = 16)


1469430 --> 1344.6
3307330 --> 50 mn ?



# I) on affecte a chaque cellule la composition du peuplement IFN
# qui a les valeurs de Dg, BA et Dprop les plus proches
#
# 1 - récupèrer Dg, BA et Dprop pour toutes les placettes LIDAR
# 2 - calculer Dg, BA et Dprop pour toutes les placettes IFN
# 3 - définir metrique de distance (travailler avec valeurs min/max des Dg, BA et Dprop?)
# 4 - mesurer les distances
# 5 - affecter la composition
#
# II) on affecte a chaque espèce un Dg et un BA qui colle avec Dg_LIDAR
# et BA_LIDAR (en utilisant les ratios de Dg et BA entre espèces observées
# sur la placette IFN)


!!!!!!! vérifier correspondance des codes TFV (factor level)
!!!!!!! pourquoi certaines cellule TFV = NA mais valeurs de Dg?





arrondir(100)
boucle tfv
créer raster ifn plot
which.min
