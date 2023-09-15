# load packages
library(tidyr)

# load loo data with assigned compo
loodf <- readRDS(paste0(tempPath, '/loodf_compo.rds'))

# load tree data and vegetation type data
if(landscape == 'bauges'){
    tree <- readRDS(paste0(tempPath, '/treeTemp.rds'))
    tree$idp <- as.character(tree$idp)
} else if (landscape == 'milicz'){
    tree <- readRDS(paste0(tempPath, '/treeTemp.rds'))
    tree$idp <- as.character(tree$idp)
    if(sum(unique(tree$idp) %in% unique(loodf$plot)) != length(unique(tree$idp))){
        stop("field and predicted data plot names do not match")
    }
} else if (landscape == 'sneznik'){
    tree <- readRDS(paste0(tempPath, '/treeTemp.rds'))
    tree$idp <- as.numeric(as.factor(tree$idp))
    tree$idp <- as.character(tree$idp)
    if(sum(unique(tree$idp) %in% unique(loodf$plot)) != length(unique(tree$idp))){
        stop("field and predicted data plot names do not match")
    }
}




###############################################################
# calculate N, Dg, BA, Dprop in NFI plots
###############################################################

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
# downscaling
###############################################################

results <- data.frame()

for (i in unique(loodf$plot)){

    # select one plot
    cell <- as.data.frame(loodf[loodf$plot == i,])

    # retrieve id composition dg and ba (LIDAR)
    id <- cell$compo
    dgtot <- cell$Dg
    batot <- cell$BA
    dprop <- cell$Dprop
    
    # deciduous - coniferous level -------------------------------------------------
    # calculate deciduous and coniferous ba (LIDAR)
    bad <- batot * dprop
    bac <- batot * (1 - dprop)

    # create dbh df
    dbh <- tree[tree$idp == id,]

    # if there is only deciduous or coniferous species in the NFI plot
    # set bad and bac to 0 and 100 accordingly
    if(nrow(dbh[dbh$spType == 'D',]) == 0){
        bad <- batot * 0
        bac <- batot * (1 - 0)
    }
    if(nrow(dbh[dbh$spType == 'C',]) == 0){
        bad <- batot * 1
        bac <- batot * (1 - 1)
    }

    # retrieve dg deciduous and coniferous (NFI)
    dgd <- as.numeric(unique(dbh[dbh$spType == 'D', 'Dgdc']))
    dgc <- as.numeric(unique(dbh[dbh$spType == 'C', 'Dgdc']))

    # calculate alpha correction coef for deciduous and coniferous
    alpha <- dgtot * sqrt( sum(bad/dgd^2, bac/dgc^2, na.rm = TRUE) / batot)

    # sp level ---------------------------------------------------------------------

    # calculate species ba (LIDAR) from NFI sp proportion
    if(landscape == 'bauges'){ # much faster on bauges landscape
        NFIplot <- NFIsp %>% filter(idp == id) %>% ungroup %>%
                            mutate(BAdclid = ifelse(spType == 'D', bad, bac),
                                    BAsplid = BAdclid * spPropdc)
    } else{
        NFIplot <- NFIsp[NFIsp$idp == id,]
        NFIplot[NFIplot$spType == 'D', 'BAdclid'] <- bad
        NFIplot[NFIplot$spType == 'C', 'BAdclid'] <- bac
        NFIplot$BAsplid <- NFIplot$BAdclid * NFIplot$spPropdc
    }

    # tree level -------------------------------------------------------------------

    # calculate tree ba (LIDAR) from NFI tree ba proportion
    dbh <- merge(dbh, NFIplot[, c('species_name', 'BAsplid')], by = 'species_name')
    dbh$BAtreelid <- dbh$BAsplid * dbh$treePropsp

    # assign a DBH (LIDAR) to all trees
    dbh$DBHlid <- dbh$DBH * alpha

    # assign weight to each tree
    dbh$wlid <- 40000 / pi * dbh$BAtreelid / dbh$DBHlid^2

    # add i and cellID now ------
    dbh$i <- i
        #   dbh$cellID25 <- cellID25
    # ---------------------------

    # calculate weight for a 25*25m pixel depending on decimal of wlid
    # could be simplified by using a Bernoulli draw (but not faster...)
    dbh$w25m <- dbh$wlid / 16
    dbh$decimal <- dbh$w25m - floor(dbh$w25m)
    dbh$random <- runif(nrow(dbh), min = 0, max = 1)
    dbh$randSmallerThanw25m <- dbh$random < dbh$decimal
    dbh[dbh$randSmallerThanw25m == TRUE, 'add'] <- 1
    dbh[dbh$randSmallerThanw25m == FALSE, 'add'] <- 0
    dbh$w25m <- floor(dbh$w25m) + dbh$add

    # remove trees with w25m = 0
    dbh <- dbh[dbh$w25m > 0,]

    # assign new dbh to trees while keeping BAtreelid
    dbh$dbhlid25m <-sqrt(40000/pi * (dbh$BAtreelid/16) / dbh$w25m)

    # return
    df <- dbh[, c('species_name', 'wlid', 'w25m', 'dbhlid25m', 'i')]
    df$wlid <- df$wlid/16
    colnames(df) <- c('sp', 'wlid', 'n', 'dbh', 'i')
    
    # if all trees are removed because of their small weights
    # an empty df must be created
    if(nrow(df) == 0){
        df <- data.frame(sp = NA, wlid = NA, n = NA, dbh = NA, i = i)
    }
    results <- rbind(results, df)

}

################################################################################
# remove trees with dbh < min dbh of tree-level inventory data
################################################################################

# define min dbh of tree-level inventory data
minDBH <- min(tree$DBH)

# remove trees < minDBH
test <- results %>% filter(dbh >= minDBH)
if(length(unique(test$i)) < length(unique(results$i))){
    eliminated <- unique(results$i)[!(unique(results$i) %in% unique(test$i))]
    results <- results %>% filter(dbh >= minDBH)
    df <- data.frame(sp = NA, wlid = NA, n = 0, dbh = 0, i = eliminated)
    results <- rbind(results, df)
}

# save
saveRDS(results, paste0(tempPath, '/loodf_results.rds'))
