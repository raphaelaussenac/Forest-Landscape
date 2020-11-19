###############################################################
# initialisation
###############################################################

# clean up environment
rm(list = ls())

# load packages
library(dplyr)
library(raster)
library(doParallel)

# set work directory
setwd("C:/Users/raphael.aussenac/Documents/GitHub/LandscapeInit")

# load sources
source('R/spTransform.R')

# load composition ID
compo <- raster("./data/Init/compoID.asc")

# load lidar data
# load LIDAR rasters
Dg <- raster("./data/Init/Dg.asc")
Dg[Dg > 1000] <- 1000 # correct the wrong min max values
BA <- raster("./data/Init/BA.asc")
BA[BA > 1000] <- 1000 # correct the wrong min max values
Dprop <- raster("./data/Init/Dprop.asc")
Dprop <- Dprop / 100

# create raster stack
rasterStack <- stack(compo, Dg, BA, Dprop)

# load NFI tree data
tree <- read.csv('./data/NFI/arbres_Bauges_2020_10_15.csv', sep = ';')

# import latin names and create deciduous / coniferous categories
tree <- spTransform(tree)

###############################################################
# calculate N, Dg, BA, Dprop in NFI plots
###############################################################

# calculate tree DBH
tree$DBH <- tree$c13 / pi
# Calculate N, Dg and BA for the whole plote, for deciduous and coniferous and
# for ech species
tree <- tree %>% mutate(BAtree =  (pi * (DBH/200)^2) * w) %>%
                 group_by(idp) %>% mutate(BAtot = sum((pi * (DBH/200)^2) * w),
                                          Dgtot = sqrt(sum(DBH^2 * w)/sum(w)),
                                          Ntot = sum(w)) %>%
                 group_by(idp, spType) %>% mutate(BAdc = sum((pi * (DBH/200)^2) * w),
                                                  Dgdc = sqrt(sum(DBH^2 * w)/sum(w)),
                                                  Ndc = sum(w)) %>%
                 group_by(idp, spType, species_name) %>% mutate(BAsp = sum((pi * (DBH/200)^2) * w),
                                                                Dgsp = sqrt(sum(DBH^2 * w)/sum(w))) %>%
                 group_by(idp, spType) %>% mutate(spPropdc = BAsp/BAdc) %>% ungroup() %>%
                 mutate(treePropsp = BAtree / BAsp)
#
# summarise at sp level
NFIsp <- tree %>% group_by(idp, species_name, BAsp, Dgsp, spType, spPropdc) %>% summarise()

###############################################################
#
###############################################################

# rasterStack <- crop(rasterStack, extent(rasterStack)/20)
# cell <- rasterStack[389]

# Split data into two separate stack rasters
rast1 <- crop(rasterStack, c(extent(rasterStack)[1],
                                  extent(rasterStack)[1] + round( (extent(rasterStack)[2] - extent(rasterStack)[1]) / 2),
                                  extent(rasterStack)[3],
                                  extent(rasterStack)[4]))
#
rast2 <- crop(rasterStack, c(extent(rasterStack)[1] + round( (extent(rasterStack)[2] - extent(rasterStack)[1]) / 2),
                                  extent(rasterStack)[2],
                                  extent(rasterStack)[3],
                                  extent(rasterStack)[4]))
#

# function to assign n trees and their dbh to each cell
assignDendro <- function(cell, i, tree, NFIsp){
  if(sum(is.na(cell)) == 0){
    # i <- 205157 # i = 205157 = plotid 4740 with 5sp = 2C + 3D:
    # deciduous - coniferous level -------------------------------------------------

    # retrieve id composition dg and ba (LIDAR)
    id <- cell[1]
    dgtot <- cell[2]
    batot <- cell[3]
    dprop <- cell[4]

    # calculate deciduous and coniferous ba (LIDAR)
    bad <- batot * dprop
    bac <- batot * (1 - dprop)

    # if there is only deciduous or coniferous species in the NFI plot
    # set bad and bac to 0 and 100 accordingly
    if(nrow(tree[tree$idp == id & tree$spType == 'D',]) == 0){
      bad <- batot * 0
      bac <- batot * (1 - 0)
    }
    if(nrow(tree[tree$idp == id & tree$spType == 'C',]) == 0){
      bad <- batot * 1
      bac <- batot * (1 - 1)
    }

    # retrieve dg deciduous and coniferous (NFI)
    dgd <- as.numeric(unique(tree[tree$idp == id & tree$spType == 'D', 'Dgdc']))
    dgc <- as.numeric(unique(tree[tree$idp == id & tree$spType == 'C', 'Dgdc']))

    # calculate alpha correction coef for deciduous and coniferous
    alphadc <- dgtot * sqrt( sum(bad/dgd^2, bac/dgc^2, na.rm = TRUE) / batot)

    # assigne corrected dg to coniferous and deciduous
    dgdcorrec <- dgd * alphadc
    dgccorrec <- dgc * alphadc

    # sp level ---------------------------------------------------------------------

    # calculate species ba (LIDAR) from NFI sp proportion
    NFIplot <- NFIsp[NFIsp$idp == id,]
    NFIplot[NFIplot$spType == 'D', 'BAdclid'] <- bad
    NFIplot[NFIplot$spType == 'C', 'BAdclid'] <- bac
    NFIplot$BAsplid <- NFIplot$BAdclid * NFIplot$spPropdc

    # calculate alpha correction coef
    NFIplot$sumRatio <- NFIplot$BAsplid / NFIplot$Dgsp^2
    # for deciduous sp
    if(nrow(NFIplot[NFIplot$spType == 'D',]) > 0){
      alphad <- dgdcorrec * sqrt( sum(NFIplot[NFIplot$spType == 'D', 'sumRatio']) / bad)
      # assign a dg (LIDAR) to all deciduous species
      NFIplot[NFIplot$spType == 'D', 'Dgsplid'] <- NFIplot[NFIplot$spType == 'D', 'Dgsp'] * alphad
    }
    # for coniferous sp
    if(nrow(NFIplot[NFIplot$spType == 'C',]) > 0){
      alphac <- dgccorrec * sqrt( sum(NFIplot[NFIplot$spType == 'C', 'sumRatio']) / bac)
      # assign a dg (LIDAR) to all coniferous species
      NFIplot[NFIplot$spType == 'C', 'Dgsplid'] <- NFIplot[NFIplot$spType == 'C', 'Dgsp'] * alphac
    }

    # calculate Nsp (LIDAR)
    NFIplot$Nsplid <- 40000/pi * ( NFIplot$BAsplid / NFIplot$Dgsplid^2 )

    # tree level -------------------------------------------------------------------

    # calculate tree ba (LIDAR) from NFI tree ba proportion
    dbh <- tree[tree$idp == id, ]
    dbh <- merge(dbh, NFIplot[, c('species_name', 'Nsplid', 'Dgsplid', 'BAsplid')], by = 'species_name')
    # dbh <- dbh %>% group_by(species_name) %>% mutate(wlid = w * Nsplid / sum(w)) %>% ungroup()
    dbh$BAtreelid <- dbh$BAsplid * dbh$treePropsp

    # calculate alpha correction coef for all trees
    dbh$sumRatio <- dbh$BAtreelid / dbh$DBH^2
    dbh$alphatree <- 999
    # dbh <- dbh %>% group_by(species_name) %>%
    #                mutate(alphatree = Dgsplid * sqrt( sum(sumRatio) / BAsplid))
    for (s in unique(dbh$species_name)){
      dbhsp <- dbh[dbh$species_name == s,]
      dbh[dbh$species_name == s, 'alphatree'] <- dbhsp$Dgsplid * sqrt( sum(dbhsp$sumRatio) / dbhsp$BAsplid)
    }
    #
    # assign a DBH (LIDAR) to all trees
    dbh$DBHlid <- dbh$DBH * dbh$alphatree

    # assign weight to each tree
    dbh$wlid <- 40000 / pi * dbh$BAtreelid / dbh$DBHlid^2

    # calculate round(weight) for a 25*25m pixel
    dbh$w25m <- max(1, round(dbh$wlid/16))

    # assign new dbh to trees while keeping BAtreelid
    dbh$dbhlid25m <-sqrt(40000/pi * (dbh$BAtreelid/16) / dbh$w25m)

    # return
    df <- dbh[, c('species_name', 'w25m', 'dbhlid25m')]
    df$i <- i
    colnames(df) <- c('sp', 'n', 'dbh', 'i')

  } else{
    df <- data.frame(sp = NA, n = NA, dbh = NA, i = i)
  }

  return(df)

}

# parallel calculation on raster cells (one raster after the other)
clustCalc <- function(rast, assignDendro, tree, NFIsp){
  # set cluster
  cl <- makeCluster(8)
  registerDoParallel(cl)
  results <- foreach(i = 1:nrow(rast[]), .combine = 'rbind', .packages = c('raster', 'rgdal')) %dopar% {assignDendro(cell = rast[i], i = i, tree, NFIsp)}
  stopCluster(cl)
  results <- results %>% arrange(i)
  return(results)
}

# run calculation
start <- Sys.time()
results <- lapply(c(rast1, rast2), clustCalc, assignDendro, tree, NFIsp)
end <- Sys.time()
end - start








# treesInit <- clustCalc(rasterStack, assignDendro, tree, NFIsp)
cl <- makeCluster(8)
registerDoParallel(cl)
results <- foreach(i = 1:nrow(rasterStack[]), .combine = 'rbind', .packages = c('raster', 'rgdal')) %dopar% {assignDendro(cell = rasterStack[i], i = i, tree, NFIsp)}
stopCluster(cl)




# verif parallel
BAinit <- results %>% group_by(i) %>% summarise(BAtot = sum((pi * (dbh/200)^2) * n * 16)) %>% arrange(i)
rasterStack$BAinit <- BAinit$BAtot
rasterStack$diff <- rasterStack$BAinit - rasterStack$BA
plot(rasterStack$diff)



plot(rasterStack$dg)
plot(rasterStack$BA)
plot(rasterStack$Dprop)
plot(rasterStack$compoID)
plot(rasterStack$BAinit)




# !!!!!!!!!!!!! difference entre arrondi et non arrondi doit etre centrée sur 0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!!!!!!!!!!!! regarder en fonction du diametre aussi !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




# verif lab !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# verifier si les BA calculer a partir des dbh lidar finaux correspondent bien
# au BAsplid, BAdclid et BAtotlid
# verifBA <- dbh %>% group_by(species_name) %>%
#                    mutate(VerifBasp = sum( (pi * (dbhlid25m/200)^2) * w25m * 16 )) %>%
#                    group_by(spType) %>%
#                    mutate(VerifBadc = sum( (pi * (dbhlid25m/200)^2) * w25m * 16 )) %>%
#                    ungroup() %>%
#                    mutate(VerifBatot = sum( (pi * (dbhlid25m/200)^2) * w25m * 16 ))
# #
# batot
# unique(verifBA$VerifBatot)
# bad
# as.numeric(unique(verifBA[verifBA$spType == 'D', 'VerifBadc']))
# bac
# as.numeric(unique(verifBA[verifBA$spType == 'C', 'VerifBadc']))
# sort(unique(verifBA$BAsplid))
# sort(unique(verifBA$VerifBasp))
#
# # calcul de Nsp a partir de Basplid et Dgsplid
# NFIplot$Nlid <- 40000/pi * ( NFIplot$BAsplid / NFIplot$Dgsplid^2 )
# # calcul de N feuillu a partir de valeurs sp
# Nspd <- sum(NFIplot[NFIplot$spType == 'D', 'Nlid'])
# Nspd
# # calcul de N coniferous a partir de valeurs sp
# Nspc <- sum(NFIplot[NFIplot$spType == 'C', 'Nlid'])
# Nspc
# # calcul de N feuillu a partir de Badlid et Dgdlid
# Ndcd <- 40000/pi * ( bad / dgdcorrec^2 )
# Ndcd
# # calcul de N coniferous a partir de Baclid et Dgclid
# Ndcc <- 40000/pi * ( bac / dgccorrec^2 )
# Ndcc
# # calcul de Ntot a partir de Batotlid et Dgtotlid
# Ntot <- 40000/pi * ( batot / dgtot^2 )
# Ntot
# Ndcd + Ndcc
# sum(NFIplot$Nlid)
# # valeurs IFN:
# NFItot[NFItot$idp == 4740,]
# verif lab !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
##################################################### other lab
# select plot with 5sp = 2C + 3D:
# nbsp <- tree %>% group_by(idp, spType) %>% summarise(n = length(unique(species_name)))
# # id 4740
# rasterTemp <- compo %in% 4740
# rasterTemp <- as.data.frame(rasterTemp)
# rasterTemp$cell <- 1:nrow(rasterTemp)
# rasterTemp <- rasterTemp[rasterTemp$layer == TRUE,]
# # e.g. 205157
# NFIplot$Dgcorrec <- NFIplot$Dgsp * alphad
# as.data.frame(NFIsp[NFIsp$idp == 4740, ])
