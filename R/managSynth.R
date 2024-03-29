managSynth <- function(landscape, sce){

  ###############################################################
  # initialisation
  ###############################################################

  # load packages
  require(dplyr)
  require(ggplot2)
  require(webr)
  require(tidyr)
  require(forcats)
  require(terra)

  # load management table depending on sce
  if(length(sce) == 1){
    df <- read.csv(paste0(landPath, '/managTableCell100_', sce, '.csv'))
  } else if(length(sce) == 2){
    df <- read.csv(paste0(landPath, '/managTableCell100_', paste0(sce[1], sce[2]), '.csv'))
  }
  df <- df %>% select(-dgModul)
  dfMaps <- df
  # load virtual tree data
  cell25 <- read.csv(paste0(landPath, '/cell25.csv'))
  # merge data and keep only forest pixels
  df <- full_join(df, cell25 %>% dplyr::select(cellID100, cellID25, forest), by = 'cellID100')
  df <- df %>% filter(forest == 1)
  # surface of a 25*25m cell in ha
  surfCell <- 0.0625

  # load cellID100 raster
  cellID100 <- rast(paste0(landPath, '/cellID100.asc'))

  ###############################################################
  # plot stand density
  ###############################################################

  if(landscape == 'bauges'){

    if(length(sce) == 1){
      pdf(paste0(evalPath, '/densityManag_', sce, '.pdf'), width = 10, height = 10)
    } else if(length(sce) == 2){
      pdf(paste0(evalPath, '/densityManag_', paste0(sce[1], sce[2]), '.pdf'), width = 10, height = 10)
    }

    pdf(paste0(evalPath, '/densityManag.pdf'), width = 10, height = 10)
    par(mfrow = c(3,2))

    # uneven conifers
    uc <- df[df$structure == 'uneven' & !is.na(df$structure) & !is.na(df$compoType) & df$compoType == 'fir and or spruce' & df$access == 1 & df$protect == 0,]
    hist(uc[uc$density == 'low', 'BA_ha'], breaks = seq(0, 110, 1), main = 'uneven conifers', xlab = 'BA_ha (red line = after logging target, dark green = abandonned)', ylab = 'nb of 25*25m cells')
    hist(uc[uc$density == 'medium', 'BA_ha'], breaks = seq(0, 110, 1), add = TRUE, col = 'green')
    hist(uc[uc$density == 'high', 'BA_ha'], breaks = seq(0, 110, 1), add = TRUE, col = 'orange')
    hist(uc[uc$manag == 'final cut', 'BA_ha'], breaks = seq(0, 110, 1), add = TRUE, col = 'green3')
    abline(v = c(25, 30, 35), lwd = 2, col = 'red')

    # uneven mixed
    comp <- c('fir and or spruce with DC', 'beech with fir and or spruce')
    um <- df[df$structure == 'uneven' & !is.na(df$structure) & !is.na(df$compoType) & df$compoType %in% comp & df$access == 1 & df$protect == 0,]
    hist(um[um$density == 'low', 'BA_ha'], breaks = seq(0, 110, 1), main = 'uneven mixed', xlab = 'BA_ha (red line = after logging target, dark green = abandonned)', ylab = 'nb of 25*25m cells')
    hist(um[um$density == 'medium', 'BA_ha'], breaks = seq(0, 110, 1), add = TRUE, col = 'green')
    hist(um[um$density == 'high', 'BA_ha'], breaks = seq(0, 110, 1), add = TRUE, col = 'orange')
    hist(um[um$manag == 'final cut', 'BA_ha'], breaks = seq(0, 110, 1), add = TRUE, col = 'green3')
    abline(v = c(20, 25, 30), lwd = 2, col = 'red')

    # uneven deciduous
    ud <- df[df$structure == 'uneven' & !is.na(df$structure) & !is.na(df$compoType) & df$compoType == 'D' & df$access == 1 & df$protect == 0,]
    hist(ud[ud$density == 'low', 'BA_ha'], breaks = seq(0, 110, 1), main = 'uneven deciduous', xlab = 'BA_ha (red line = after logging target, dark green = abandonned)', ylab = 'nb of 25*25m cells')
    hist(ud[ud$density == 'medium', 'BA_ha'], breaks = seq(0, 110, 1), add = TRUE, col = 'green')
    hist(ud[ud$density == 'high', 'BA_ha'], breaks = seq(0, 110, 1), add = TRUE, col = 'orange')
    hist(ud[ud$manag == 'final cut', 'BA_ha'], breaks = seq(0, 110, 1), add = TRUE, col = 'green3')
    abline(v = c(15, 20, 25), lwd = 2, col = 'red')

    # even except deciduous
    even <- df[df$structure == 'even' & !is.na(df$structure) & df$access == 1 & df$protect == 0 & df$compoType != 'D' & !is.na(df$compoType), ]
    hist(even[even$density == 'medium', 'rdi'], breaks = seq(0, 2, 0.05), main = 'even except deciduous', xlab = 'rdi (red line = after logging target, dark green = abandonned)', ylab = 'nb of 25*25m cells')
    hist(even[even$density == 'high', 'rdi'], breaks = seq(0, 2, 0.05), add = TRUE, col = 'orange')
    hist(even[even$manag == 'final cut', 'rdi'], breaks = seq(0, 2, 0.05), add = TRUE, col = 'green3')
    abline(v = c(0.55, 0.65), lwd = 2, col = 'red')

    # even deciduous
    evenD <- df[df$structure == 'even' & !is.na(df$structure) & df$access == 1 & df$protect == 0 & df$compoType == 'D' & !is.na(df$compoType), ]
    plot(evenD$Dg, evenD$rdi, pch = 16, main = 'even deciduous', xlab = 'Dg (red = coppice, dark green = abandonned)', ylab = 'rdi')
    points(evenD[evenD$manag == 'coppice', 'Dg'], evenD[evenD$manag == 'coppice', 'rdi'], col = 'red', pch = 16)
    points(evenD[evenD$manag == 'final cut', 'Dg'], evenD[evenD$manag == 'final cut', 'rdi'], col = 'green3', pch = 16)

    par(mfrow = c(1,1))
    dev.off()

  } else if(landscape == 'milicz' | landscape == 'sneznik'){

    # BA_ha distribution
    plotBA <- ggplot(data = df %>% filter(!is.na(compoType))) +
    geom_histogram(aes(BA_ha), bins = 100) +
    facet_wrap(.~compoType) +
    theme_bw() +
    ylab('nb of 25*25m cells')
    plotBA
    if(length(sce) == 1){
      ggsave(file = paste0(evalPath, '/densityBA_', sce, '.pdf'), plot = plotBA, width = 10, height = 10)
    } else if(length(sce) == 2){
      ggsave(file = paste0(evalPath, '/densityBA_', paste0(sce[1], sce[2]), '.pdf'), plot = plotBA, width = 10, height = 10)
    }

    # RDI distribution
    plotRDI <- ggplot(data = df %>% filter(!is.na(compoType))) +
    geom_histogram(aes(rdi), bins = 100) +
    facet_wrap(.~compoType) +
    theme_bw() +
    ylab('nb of 25*25m cells')
    plotRDI
    if(length(sce) == 1){
      ggsave(file = paste0(evalPath, '/densityRDI_', sce, '.pdf'), plot = plotRDI, width = 10, height = 10)
    } else if(length(sce) == 2){
      ggsave(file = paste0(evalPath, '/densityRDI_', paste0(sce[1], sce[2]), '.pdf'), plot = plotRDI, width = 10, height = 10)
    }
  }


  ###############################################################
  # plot composition types proportion in landscape
  ###############################################################

  compo <- df %>% filter(!is.na(compoType)) %>% group_by(compoType) %>% mutate(surface = surfCell) %>%
                      summarise(surface = sum(surface)) %>% arrange(-surface) %>%
                      ungroup() %>% mutate(totSurf = sum(surface)) %>%
                      group_by(compoType) %>% mutate(relSurf = surface * 100 / totSurf) %>%
                      dplyr::select(-totSurf)
  compo$compoType <- factor(compo$compoType, levels = compo$compoType)

  # plot composition at 1ha cells
  pl1 <- ggplot(data = compo) +
  geom_bar(aes(x = compoType, y = surface), stat = 'identity') +
  geom_text(aes(x = compoType, y = surface, label = paste(round(relSurf, 2), '%') ), vjust = -0.3, size = 5, col = 'orange') +
  geom_text(aes(x = compoType, y = surface, label = paste(round(surface, 2), 'ha') ), vjust = +1.2, size = 5, col = 'orange') +
  theme_light() +
  xlab('composition') +
  ylab('surface (ha)') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
       plot.title = element_text(hjust = 0.5))
  pl1
  if(length(sce) == 1){
    ggsave(file = paste0(evalPath, '/compoType_', sce, '.pdf'), plot = pl1, width = 10, height = 10)
  } else if(length(sce) == 2){
    ggsave(file = paste0(evalPath, '/compoType_', paste0(sce[1], sce[2]), '.pdf'), plot = pl1, width = 10, height = 10)
  }

  ###############################################################
  # synthesise management type for each composition
  ###############################################################

  if(landscape == 'bauges'){

    # funcrtion to make sure all tables have similar structure
    # in order to rbind them
    fillTab <- function(tab){
      # create table structure
      structureTab <- data.frame(matrix(ncol = 9, nrow = nrow(tab)))
      names(structureTab) <- c('access', 'owner', 'structure', 'density', 'beech with fir and or spruce', 'D', 'fir and or spruce', 'fir and or spruce with DC', 'other compo')
      # fill the table
      for (i in names(tab)){
        structureTab[, i] <- tab[, i]
      }
      return(structureTab)
    }

    # accessible stands
    standType <- df
    standType[standType$protect == 1, 'access'] <- 0
      # [e]ven and [u]neven [s]tands
    eus <- standType %>% filter(!is.na(compoType), !(manag %in% c('coppice', 'final cut'))) %>% group_by(access, compoType, owner, structure, density) %>% summarise(surf = n()*surfCell) %>% ungroup()
    tab2 <- eus %>% filter(access == 1) %>% pivot_wider(names_from = compoType, values_from = surf)
    tab2 <- fillTab(tab2)
     # [cop]pice stands
    cop <- standType %>% filter(!is.na(compoType), manag == 'coppice') %>% group_by(access, compoType, owner, structure) %>% summarise(surf = n()*surfCell) %>% ungroup()
    tab3 <- cop %>% filter(access == 1) %>% pivot_wider(names_from = compoType, values_from = surf)
    tab3$structure <- 'coppice'
    tab3$density <- NA
    tab3 <- fillTab(tab3)
     # [fin]al cut stands
    fin <- standType %>% filter(!is.na(compoType), manag == 'final cut') %>% group_by(access, compoType, owner) %>% summarise(surf = n()*surfCell) %>% ungroup()
    tab4 <- fin %>% filter(access == 1) %>% pivot_wider(names_from = compoType, values_from = surf)
    tab4$structure <- 'final cut'
    tab4$density <- NA
    tab4 <- fillTab(tab4)
    # [inacc]essible stands
    inacc <- standType %>% filter(!is.na(compoType)) %>% group_by(access, compoType, owner, structure, density) %>% summarise(surf = n()*surfCell) %>% ungroup()
    tab1 <- inacc %>% filter(access == 0) %>% group_by(compoType, owner) %>%
                      summarise(surf = sum(surf)) %>%
                      pivot_wider(names_from = compoType, values_from = surf) %>%
                      mutate(access = 0, structure = NA, density = NA)
    tab1 <- fillTab(tab1)
    # rbind tables
    tab <- rbind(tab1, tab3, tab4, tab2)
    tab <- tab %>% rename(accessProtec = access, manag = structure)

  } else if(landscape == 'milicz'){
    tab <- df %>% group_by(manag, compoType) %>% summarise(surf = n()*surfCell) %>%
                  filter(!is.na(compoType))
    tab1 <- tab %>% filter(manag == 'no manag') %>% pivot_wider(names_from = compoType, values_from = surf)
    missingCol <- unique(tab$compoType)[!(unique(tab$compoType) %in% names(tab1))]
    tab1[, missingCol] <-  0
    tab2 <- tab %>% filter(manag != 'no manag') %>% pivot_wider(names_from = compoType, values_from = surf)
    missingCol <- unique(tab$compoType)[!(unique(tab$compoType) %in% names(tab2))]
    tab2[, missingCol] <-  0
    tab <- rbind(tab1[, names(tab2)], tab2) %>% ungroup()
  } else if(landscape == 'sneznik'){
    tab <- df %>% group_by(manag, compoType) %>% summarise(surf = n()*surfCell) %>% pivot_wider(names_from = compoType, values_from = surf)
  }

  if(length(sce) == 1){
    write.csv(tab, paste0(evalPath, '/managSyn_', sce, '.csv'), row.names = F)
  } else if(length(sce) == 2){
    write.csv(tab, paste0(evalPath, '/managSyn_', paste0(sce[1], sce[2]), '.csv'), row.names = F)
  }


  ###############################################################
  # surface barplot of each type
  ###############################################################

  if(landscape == 'bauges'){
    stand <- df %>% filter(!is.na(forestCellsPerHa))
    stand[stand$manag %in% c('coppice', 'final cut'), 'structure'] <- stand[stand$manag %in% c('coppice', 'final cut'), 'manag']
    stand[stand$access == 0 | stand$protect == 1, 'structure'] <- 'inacc&protec'
    stand <- stand %>% group_by(compoType, structure) %>% summarise(surf = n()*surfCell) %>% arrange(-surf)
    stand$structure <- factor(stand$structure, levels = c('uneven', 'inacc&protec', 'even', 'final cut', 'coppice'))
    stand$compoType <- factor(stand$compoType, levels = c('fir and or spruce with DC', 'D with fir and or spruce', 'beech with fir and or spruce', 'fir and or spruce', 'D', 'DC with fir and or spruce', 'beech', 'C with fir and or spruce', 'DC', 'other compo'))
  } else if(landscape == 'milicz' | landscape == 'sneznik'){
    stand <- df %>% filter(!is.na(forestCellsPerHa))
    stand[stand$access == 0 | stand$protect == 1, 'structure'] <- 'inacc&protec'
    stand <- stand %>% group_by(compoType, structure) %>% summarise(surf = n()*surfCell) %>% arrange(-surf)
    stand$structure <- factor(stand$structure, levels = c('even', 'uneven', 'inacc&protec'))
    stand$compoType <- fct_reorder(stand$compoType, stand$surf, max, .desc = TRUE)
  }

  # plot
  pl2 <- ggplot() +
  geom_bar(data = stand, aes(x = compoType, y = surf, fill = structure), stat = 'identity', position = position_dodge(width = 1)) +
  # scale_fill_manual(values = c('green4', 'orange4', 'green3', 'blue1', 'orange')) +
  theme_minimal()
  pl2
  if(length(sce) == 1){
    ggsave(file = paste0(evalPath, '/managSurf1_', sce, '.pdf'), plot = pl2, width = 15, height = 10)
  } else if(length(sce) == 2){
    ggsave(file = paste0(evalPath, '/managSurf1_', paste0(sce[1], sce[2]), '.pdf'), plot = pl2, width = 15, height = 10)
  }

  ###############################################################
  # piedonut of surface
  ###############################################################

  if(length(sce) == 1){
    pdf(paste0(evalPath, '/managSurf2_', sce, '.pdf'), width = 10, height = 10)
  } else if(length(sce) == 2){
    pdf(paste0(evalPath, '/managSurf2_', paste0(sce[1], sce[2]), '.pdf'), width = 10, height = 10)
  }

  if(landscape == 'bauges'){
    stand1 <- df
    stand1[stand1$protect == 1 & !is.na(stand1$protect), 'access'] <- 0
    stand1[stand1$manag %in% c('coppice', 'final cut'), 'structure'] <- stand1[stand1$manag %in% c('coppice', 'final cut'), 'manag']
    stand1 <- stand1 %>% filter(!is.na(compoType)) %>% group_by(access, compoType, owner, structure, density) %>% summarise(surf = n()*surfCell) %>% ungroup()
    stand1$compoType <- factor(stand1$compoType, levels = c('fir and or spruce with DC', 'D with fir and or spruce', 'beech with fir and or spruce', 'fir and or spruce', 'D', 'DC with fir and or spruce', 'beech', 'C with fir and or spruce', 'DC', 'other compo'))
    stand1$compoType <- droplevels(stand1$compoType)
    stand1$structure <- as.character(stand1$structure)
    stand1[stand1$access == 0, 'structure'] <- 'inacc/protect'
    stand1$structure <- factor(stand1$structure, levels = c('inacc/protect', 'uneven', 'even', 'final cut', 'coppice'))
    donut <- stand1 %>% group_by(compoType, structure) %>% summarise(surf = sum(surf)) %>% ungroup() %>%
                      arrange(compoType, structure, -surf) %>% mutate(ymax = cumsum(surf),
                                                ymin = lag(ymax, default = 0))
    #
    PieDonut(donut, aes(compoType, structure, count = surf), showPieName = FALSE, start = pi/2, title = 'inacc/protect / uneven / even / final cut / coppice')
    dev.off()

  } else if(landscape %in% c('milicz', 'sneznik')){

    stand1 <- df
    stand1 <- stand1 %>% filter(!is.na(compoType)) %>% group_by(compoType, manag) %>% summarise(surf = n()*surfCell) %>% ungroup()
    stand1$compoType <- fct_reorder(stand1$compoType, stand1$surf, max, .desc = TRUE)
    stand1$manag <- fct_reorder(stand1$manag, stand1$surf, max, .desc = TRUE)
    donut <- stand1 %>% group_by(compoType, manag) %>% summarise(surf = sum(surf)) %>% ungroup() %>%
                      arrange(compoType, manag, -surf) %>% mutate(ymax = cumsum(surf),
                                                ymin = lag(ymax, default = 0))
    #
    PieDonut(donut, aes(compoType, manag, count = surf), showPieName = FALSE, start = pi/2, title = 'even - uneven / protect OR no manag')
    dev.off()

  }


  ###############################################################
  # maps
  ###############################################################

  # management map
  # add manag to cellID100 raster
  dfMaps$manag <- as.factor(dfMaps$manag)
  cellID100$manag <- as.numeric(dfMaps$manag)
  plot(cellID100$manag)
  # convert raster into polygon
  managPoly <- terra::as.polygons(cellID100$manag, trunc = TRUE, dissolve = TRUE, na.rm = TRUE)
  # assign management name to polygons (insted of intergers codes)
  corr <- data.frame(level = c(1:length(levels(dfMaps$manag))), levels(dfMaps$manag))
  managPoly <- merge(managPoly, corr, by.x = 'manag', by.y = 'level')
  managPoly$manag <- NULL
  names(managPoly) <- 'manag'
  # save
  if(length(sce) == 1){
    writeVector(managPoly, paste0(evalPath, '/managPoly_', sce, '.shp'), overwrite=TRUE)
  } else if(length(sce) == 2){
    writeVector(managPoly, paste0(evalPath, '/managPoly_', paste0(sce[1], sce[2]), '.shp'), overwrite=TRUE)
  }




  # composition type map
  dfMaps$compoType <- as.factor(dfMaps$compoType)
  cellID100$compoType <- as.numeric(dfMaps$compoType)
  plot(cellID100$compoType)
  # convert raster into polygon
  compoTypePoly <- terra::as.polygons(cellID100$compoType, trunc = TRUE, dissolve = TRUE, na.rm = TRUE)
  # assign management name to polygons (insted of intergers codes)
  corr <- data.frame(level = c(1:length(levels(dfMaps$compoType))), levels(dfMaps$compoType))
  compoTypePoly <- merge(compoTypePoly, corr, by.x = 'compoType', by.y = 'level')
  compoTypePoly$compoType <- NULL
  names(compoTypePoly) <- 'compoType'
  # save

  if(length(sce) == 1){
    writeVector(compoTypePoly, paste0(evalPath, '/compoTypePoly_', sce, '.shp'), overwrite=TRUE)
  } else if(length(sce) == 2){
    writeVector(compoTypePoly, paste0(evalPath, '/compoTypePoly_', paste0(sce[1], sce[2]), '.shp'), overwrite=TRUE)
  }

}


# create synthesis for all alternative managements
altManagSynth <- function(landscape){
  lapply(c('B', 'I', 'E'), managSynth, landscape = landscape)
  managSynth(landscape, sce = c('B', 'C'))
  managSynth(landscape, sce = c('I', 'C'))
  managSynth(landscape, sce = c('E', 'C'))
}
# sce = c('I', 'E', 'B', 'C')
# for [I]ntensification / [E]xtensification / [B]aseline / working for [C]omplexity
