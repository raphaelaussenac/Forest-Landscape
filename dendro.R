###############################################################
# initialisation
###############################################################

# clean up environment
rm(list = ls())

# load packages
library(raster)
library(dplyr)

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
# calculate Dg, BA, Dprop... for NFI plots
###############################################################

# calculate tree DBH
tree$DBH <- tree$c13 / pi
#
# calculate total Dg and BA
NFItot <- tree %>% group_by(idp) %>% summarise(BAtot = sum((pi * (DBH/200)^2) * w),
                                               Dgtot = sqrt(sum(DBH^2 * w)/sum(w)),
                                               N = sum(w))
#
# calculate Dg and BA for deciduous and coniferous
NFIdc <- tree %>% group_by(idp, spType) %>% summarise(BAdc = sum((pi * (DBH/200)^2) * w),
                                                      Dgdc = sqrt(sum(DBH^2 * w)/sum(w)),
                                                      N = sum(w))
#
# calculate Dg and BA for all species
NFIsp <- tree %>% group_by(idp, spType, species_name) %>% summarise(BAsp = sum((pi * (DBH/200)^2) * w),
                                                      Dgsp = sqrt(sum(DBH^2 * w)/sum(w))) %>%
                  group_by(idp, spType) %>% mutate(BAdc = sum(BAsp)) %>%
                  mutate(spPropdc = BAsp/BAdc)
#

###############################################################
#
###############################################################

i <- 205157
# deciduous - coniferous level -------------------------------------------------

# retrieve id composition dg and ba (LIDAR)
id <- rasterStack[i][1]
dgtot <- rasterStack[i][2]
batot <- rasterStack[i][3]
dprop <- rasterStack[i][4]

# calculate deciduous and coniferous ba (LIDAR)
bad <- batot * dprop
bac <- batot * (1 - dprop)

# retrieve dg deciduous and coniferous (NFI)
dgd <- as.numeric(NFIdc[NFIdc$idp == id & NFIdc$spType == 'D', 'Dgdc'])
dgc <- as.numeric(NFIdc[NFIdc$idp == id & NFIdc$spType == 'C', 'Dgdc'])

# calculate alpha correction coef for deciduous and coniferous
alphadc <- dgtot * sqrt( sum(bad/dgd^2, bac/dgc^2) / batot)

# assigne corrected dg to coniferous and deciduous
dgdcorrec <- dgd * alphadc
dgccorrec <- dgc * alphadc


# sp level ---------------------------------------------------------------------

# calculate species ba (LIDAR) from NFI sp proportion
NFIplot <- NFIsp[NFIsp$idp == id,]
NFIplot[NFIplot$spType == 'D', 'BAdclid'] <- bad
NFIplot[NFIplot$spType == 'C', 'BAdclid'] <- bac
NFIplot$BAsplid <- NFIplot$BAdclid * NFIplot$spPropdc
NFIplot$sumRatio <- NFIplot$BAsplid / NFIplot$Dgsp^2

# calculate alpha correction coef for deciduous sp
alphad <- dgdcorrec * sqrt( sum(NFIplot[NFIplot$spType == 'D', 'sumRatio']) / bad)
# calculate alpha correction coef for coniferous sp
alphac <- dgccorrec * sqrt( sum(NFIplot[NFIplot$spType == 'C', 'sumRatio']) / bac)

# assign a dglid to all species
NFIplot$Dgsplid <- NFIplot$Dgsp * alphad


# verif lab !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

# III) alpha pour les diamètres ------------------------------------------------




###################################################### lab

nbsp <- tree %>% group_by(idp, spType) %>% summarise(n = length(unique(species_name)))
# 5sp = 2C + 3D:
# id 4740
rasterTemp <- compo %in% 4740
rasterTemp <- as.data.frame(rasterTemp)
rasterTemp$cell <- 1:nrow(rasterTemp)
rasterTemp <- rasterTemp[rasterTemp$layer == TRUE,]
# e.g. 205157
NFIplot$Dgcorrec <- NFIplot$Dgsp * alphad


as.data.frame(NFIdc[NFIdc$idp == 4740, ])
