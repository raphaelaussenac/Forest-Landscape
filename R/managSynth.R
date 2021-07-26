managSynth <- function(){

  ###############################################################
  # initialisation
  ###############################################################

  # load packages
  require(dplyr)
  require(ggplot2)
  library(webr)
  library(tidyr)

  # load management table
  df <- read.csv(paste0(landPath, '/managTable.csv'))

  ###############################################################
  # plot composition types proportion in landscape
  ###############################################################

  compo <- df %>% filter(!is.na(compoType)) %>% group_by(compoType) %>% mutate(surface = 1) %>%
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
  xlab('main species') +
  ylab('surface (ha)') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
       plot.title = element_text(hjust = 0.5))
  pl1

  ggsave(file = paste0(evalPath, '/compoType.pdf'), plot = pl1, width = 10, height = 10)

  ###############################################################
  # synthesise management type for each composition
  ###############################################################

  standType <- df
  standType[standType$protect == 1, 'access'] <- 0
  standType <- standType %>% filter(!is.na(compoType)) %>% group_by(access, compoType, owner, structure, density) %>% summarise(surf = n()) %>% ungroup()
  standType$structure <- factor(standType$structure, levels = c('uneven', 'even'))
  standType$compoType <- factor(standType$compoType, levels = c('D with fir and/or spruce', 'beech with fir and/or spruce', 'fir and/or spruce', 'D', 'DC with fir and/or spruce', 'beech', 'C with fir and/or spruce', 'DC'))
  tab2 <- standType %>% filter(access == 1) %>% pivot_wider(names_from = compoType, values_from = surf)
  tab1 <- standType %>% filter(access == 0) %>% group_by(compoType, owner) %>%
                    summarise(surf = sum(surf)) %>%
                    pivot_wider(names_from = compoType, values_from = surf) %>%
                    mutate(access = 0, structure = NA, density = NA)
  tab1 <- tab1[, names(tab2)]
  tab <- rbind(tab1, tab2)
  tab <- tab %>% rename(accessProtec = access)

  write.csv(tab, paste0(evalPath, '/managSyn.csv'), row.names = F)

  ###############################################################
  # surface barplot of each type
  ###############################################################

  stand <- df %>% filter(!is.na(compoType)) %>% group_by(protect, access, compoType, owner, structure, density) %>% summarise(surf = n()) %>% ungroup()
  stand$structure <- factor(stand$structure, levels = c('uneven', 'even', 'inaccessible', 'protected'))
  stand$compoType <- factor(stand$compoType, levels = c('D with fir and/or spruce', 'beech with fir and/or spruce', 'fir and/or spruce', 'D', 'DC with fir and/or spruce', 'beech', 'C with fir and/or spruce', 'DC'))
  stand$structure2 <- stand$structure
  mainType <- stand %>% group_by(compoType, structure) %>% summarise(surf = sum(surf)) %>% arrange(-surf)
  mainType$structure2 <- mainType$structure
  noAcc <- stand %>% filter(access == 0) %>% group_by(compoType, structure) %>% summarise(surf = sum(surf))
  noAcc$structure2 <- 'inaccessible'
  protect <- stand %>% filter(protect == 1) %>% group_by(compoType, structure) %>% summarise(surf = sum(surf))
  protect$structure2 <- 'protected'

  # plot
  pl2 <- ggplot() +
  geom_bar(data = mainType, aes(x = compoType, y = surf, group = structure, fill = structure2), stat = 'identity', position = position_dodge(width = 1)) +
  geom_bar(data = noAcc, aes(x = compoType, y = surf, group = structure, fill = structure2), width = 0.5, stat = 'identity', position = position_dodge(width = 1)) +
  geom_bar(data = protect, aes(x = compoType, y = surf, group = structure, fill = structure2), width = 0.25, stat = 'identity', position = position_dodge(width = 1)) +
  scale_fill_manual(values = c("green3", "orange3", "red", 'green4')) +
  theme_minimal()
  pl2

  ggsave(file = paste0(evalPath, '/managSurf1.pdf'), plot = pl2, width = 15, height = 10)

  ###############################################################
  # piedonut of surface
  ###############################################################

  stand1 <- df
  stand1[stand1$protect == 1 & !is.na(stand1$protect), 'access'] <- 0
  stand1 <- stand1 %>% filter(!is.na(compoType)) %>% group_by(access, compoType, owner, structure, density) %>% summarise(surf = n()) %>% ungroup()
  stand1$compoType <- factor(stand1$compoType, levels = c('D with fir and/or spruce', 'beech with fir and/or spruce', 'fir and/or spruce', 'D', 'DC with fir and/or spruce', 'beech', 'C with fir and/or spruce', 'DC'))
  stand1$structure <- as.character(stand1$structure)
  stand1[stand1$access == 0, 'structure'] <- 'inacc/protect'
  stand1$structure <- factor(stand1$structure, levels = c('inacc/protect', 'uneven', 'even'))
  donut <- stand1 %>% group_by(compoType, structure) %>% summarise(surf = sum(surf)) %>% ungroup() %>%
                    arrange(compoType, structure, -surf) %>% mutate(ymax = cumsum(surf),
                                              ymin = lag(ymax, default = 0))
  #
  pdf(paste0(evalPath, '/managSurf2.pdf'), width = 10, height = 10)
  PieDonut(donut, aes(compoType, structure, count = surf), showPieName = FALSE, start = pi/2)
  dev.off()


  ###############################################################
  # management map
  ###############################################################

  # cellID100$standType <- as.numeric(as.factor(df$standType))
  # plot(cellID100$standType)
  # plot(park, add = T)

  # TODO: verifier les cartes

}
