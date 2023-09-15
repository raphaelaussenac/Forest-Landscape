# load packages
require(stringr)
require(dplyr)



# load loo data
load('./loo/LOO_values.RDA')
if(landscape == 'bauges'){
    loodf <- Bauges
    load('./loo/bauges_tfv.RDA')
    tfv <- Bauges %>% select(plotMerge, CODE_TFV)
    loodf <- left_join(loodf, tfv, join_by('plot' == 'plotMerge'))
    # load tree data and vegetation type data
    tree <- readRDS(paste0(tempPath, '/treeTemp.rds'))
    tree$idp <- as.character(tree$idp)
} else if (landscape == 'milicz'){
    loodf <- Milicz %>% mutate(CODE_TFV = 1)
    # load tree data and vegetation type data
    tree <- readRDS(paste0(tempPath, '/treeTemp.rds'))
    tree$idp <- as.character(tree$idp)
    if(sum(unique(tree$idp) %in% unique(loodf$plot)) != length(unique(tree$idp))){
        stop("field and predicted data plot names do not match")
    }
} else if (landscape == 'sneznik'){
    loodf <- Sneznik %>% mutate(CODE_TFV = 1)
    loodf$plot <- as.character(as.numeric(as.factor(loodf$plot)))
    # load tree data and vegetation type data
    tree <- readRDS(paste0(tempPath, '/treeTemp.rds'))
    tree$idp <- as.character(tree$idp)
    if(sum(unique(tree$idp) %in% unique(loodf$plot)) != length(unique(tree$idp))){
        stop("field and predicted data plot names do not match")
    }
}
loodf <- loodf %>% select(plot, names(loodf)[str_detect(names(loodf), 'predicted')], CODE_TFV) %>%
                   rename(BA = BA_predicted, Dg = Dg_predicted, Dprop = DP_predicted)





###############################################################
# group TFV types together
###############################################################

if(landscape == 'bauges'){
    # main types
    A <- 'FF1-00-00' #: Forêt fermée à mélange de feuillus
    B <- 'FF2G61-61' #: Forêt fermée de sapin ou épicéa
    C <- 'FF31' #: Forêt fermée à mélange de feuillus prépondérants et conifères
    D <- 'FF32' #: Forêt fermée à mélange de conifères prépondérants et feuillus
    E <- 'FF1-09-09' #: Forêt fermée de hêtre pur

    groupTfv <- function(df){
        df <- df[!(df$CODE_TFV %in% c('FO0', 'FF0', 'LA4', 'LA6')) & !is.na(df$CODE_TFV),]
        df[df$CODE_TFV == 'FO1', 'CODE_TFV'] <- A
        df[df$CODE_TFV == 'FF1-00', 'CODE_TFV'] <- A
        df[df$CODE_TFV == 'FO2', 'CODE_TFV'] <- B
        df[df$CODE_TFV == 'FO3', 'CODE_TFV'] <- C
        df[df$CODE_TFV == 'FF1-49-49', 'CODE_TFV'] <- A
        df[df$CODE_TFV == 'FF2-00-00', 'CODE_TFV'] <- B
        df[df$CODE_TFV == 'FF1G01-01', 'CODE_TFV'] <- A
        df[df$CODE_TFV == 'FF2G53-53', 'CODE_TFV'] <- B
        df[df$CODE_TFV == 'FF2-90-90', 'CODE_TFV'] <- B
        df[df$CODE_TFV == 'FF2-64-64', 'CODE_TFV'] <- B
        df[df$CODE_TFV == 'FF2-00', 'CODE_TFV'] <- B
        df[df$CODE_TFV == 'FF1-10-10', 'CODE_TFV'] <- A
        df[df$CODE_TFV == 'FF2-52-52', 'CODE_TFV'] <- B
        df[df$CODE_TFV == '', 'CODE_TFV'] <- A
        return(df)
    }

# assign new TFV to spatial polygons
loodf <- groupTfv(loodf)
tree <- groupTfv(tree)

}




###############################################################
# Get Dg, BA, Dprop and TFV of all calibration plots
###############################################################

NFIDprop <- tree %>% group_by(idp, spType) %>%
                            summarise(BAdc = sum((pi * (DBH/200)^2) * w)) %>%
                            group_by(idp) %>% mutate(BA = sum(BAdc), Dprop = BAdc/BA) %>%
                            filter(spType == 'D') %>% select(-spType, -BAdc, -BA)
# calculate Dg and BA
NFI <- tree %>% group_by(idp) %>%
                      summarise(Dg = sqrt(sum(DBH^2 * w)/sum(w)),
                                BA = sum((pi * (DBH/200)^2) * w),
                                CODE_TFV = unique(CODE_TFV))
NFI <- merge(NFIDprop, NFI, by = 'idp', all = TRUE)
NFI[is.na(NFI$Dprop), 'Dprop'] <- 0

# convert NFI TFV CODE in numeric values
if (landscape == 'bauges'){
    NFI$CODE_TFV <- factor(NFI$CODE_TFV, levels = c('FF1-00-00', 'FF2G61-61', 'FF31', 'FF32', 'FF1-09-09'))
    loodf$CODE_TFV <- factor(loodf$CODE_TFV, levels = c('FF1-00-00', 'FF2G61-61', 'FF31', 'FF32', 'FF1-09-09'))
}
NFI$CODE_TFV <- as.numeric(NFI$CODE_TFV)
loodf$CODE_TFV <- as.numeric(loodf$CODE_TFV)


# create column for nearest plot id = compoID
loodf$compo <- NA

# Assign compo to each plot
for (i in 1:nrow(loodf)){
    
    # select one plot
    cell <- as.data.frame(loodf[i,])

    # exclude (loo) NFI plots in Sneznik and Milicz
    pl <- cell$plot
    NFI2 <- NFI %>% filter(idp != pl)

    # normalise Dg and BA values to 0-1 range
    minDg <- min(min(NFI2$Dg), min(cell$Dg))
    maxDg <- max(max(NFI2$Dg), max(cell$Dg))
    NFI2$Dg01 <- ( NFI2$Dg - minDg ) / ( maxDg - minDg )
    cell$Dg01 <- ( cell$Dg - minDg ) / ( maxDg - minDg )
    minBA <- min(min(NFI2$BA), min(cell$BA))
    maxBA <- max(max(NFI2$BA), max(cell$BA))
    NFI2$BA01 <- ( NFI2$BA - minBA ) / ( maxBA - minBA )
    cell$BA01 <- ( cell$BA - minBA ) / ( maxBA - minBA )
    NFI2 <- NFI2 %>% select(idp, Dprop, Dg01, BA01, CODE_TFV)

    # get list of plots (and their Dg, BA, Dprop values) in the TFV
    # type associated to the cell
    plotlist <- NFI2[NFI2$CODE_TFV == cell[, 'CODE_TFV'], ]
    # get cell values = coordinates
    plotlist$DpropCell <- cell[, 'Dprop']
    plotlist$Dg01Cell <- cell[, 'Dg01']
    plotlist$BA01Cell <- cell[, 'BA01'] 
    # calculate distance
    plotlist$distance <- sqrt( (plotlist$Dprop - plotlist$DpropCell)^2 + (plotlist$Dg01 - plotlist$Dg01Cell)^2 + (plotlist$BA01 - plotlist$BA01Cell)^2)
    # determine nearest NFI plot
    NearestPlot <- plotlist[plotlist$distance == min(plotlist$distance), 'idp']
    NearestPlot <- NearestPlot[1]
    # save in loodf
    loodf[loodf$plot == cell$plot, 'compo'] <- NearestPlot
}

# print error message if nearest plot missing
if(sum(is.na(loodf$compo)) > 0){
  stop("nearest plot missing")
}

# save
saveRDS(loodf, paste0(tempPath, '/loodf_compo.rds'))
