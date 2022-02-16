evalHeight <- function(landscape){

  ###############################################################
  # initialisation
  ###############################################################

  # load packages
  require(raster)
  require(dplyr)
  require(ggplot2)

  # load LiDAR heights
  lidH <- list.files(paste0('./data/', landscape, '/LiDARh'), pattern = '.rda')
  load(paste0('./data/', landscape, '/LiDARh/', lidH))

  # load tree data
  tree <- read.csv(paste0(landPath, '/trees.csv'))

  # load cellID25 raster
  cellID25 <- raster(paste0(landPath, '/cellID25.asc'))

  # merge data sets
  # stack cellID25 and LiDAR heights
  rastH <- stack(cellID25, metrics.map.Hm)
  # convert into data.frame
  lidhdf <- as.data.frame(rastH)
  # merge with tree data
  tree <- left_join(tree, lidhdf, by = 'cellID25')

  ###############################################################
  # compare tree heights and LiDAR heights
  ###############################################################

  # calculate mean height of n highest trees
  compareH <- function(nbTrees, tree){
    compH <- tree %>% group_by(cellID25) %>% arrange(-cellID25, -h) %>% slice(1:nbTrees) %>%
                      mutate(nh = n()) %>% filter(N>=nbTrees, nh == nbTrees) %>%
                      summarise(h = mean(h), hlid = unique(get(paste0('Hm', nbTrees))))
  }
  heights <- compareH(6, tree)
  mod <- lm(hlid~h, data = heights)
  summary(mod)

  # create classes for h --> bowplot
  heights <- heights %>% mutate(bin=cut_width(h, width = 1, center = 0)) %>%
                         group_by(bin) %>% mutate(class = round(mean(h))) %>% ungroup()
  # heights <- heights %>% mutate(bin=cut(h, 40, labels = FALSE))

  # plot
  pl1 <- ggplot(heights, aes(x = h, y = hlid)) +
  geom_point(alpha = 0.5, size = 0.5, col = 'grey', pch = 16) +
  xlim(5, 45) +
  ylim(5, 45) +
  coord_fixed() +
  annotate(geom = 'text', x = 10, y = 40, label = paste('r² = ', round(summary(mod)$r.squared,2)), col = 'red', size = 10) +
  ggtitle(landscape) +
  geom_boxplot(aes(x = class, y = hlid, group = class), alpha = 0.8) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = 2, size = 1) +
  geom_abline(intercept = mod$coef[1], slope = mod$coef[2], color = "red", linetype = 1, size = 1) +
  theme_bw()
  pl1
  ggsave(file = paste0(evalHeightPath, '/', landscape, '_evalHeight.jpg'), plot = pl1, width = 10, height = 10)

  # TODO: JMM a mis des limites min de diametre ????
  # TODO: inclure toutes les cellules = comparaison de hdom de max 6 arbres
  # TODO: ajouter nb/proportion de cellules sur le territoire (ou hlid & h > 6)

  ###############################################################
  # mean square deviation
  ###############################################################

  # Y (observations)
  colnames(heights)[colnames(heights) == 'hlid'] <- 'Y'
  # X (predictions)
  colnames(heights)[colnames(heights) == 'h'] <- 'X'

  # MSD calculation
  heights <- heights %>% dplyr::select(-bin, -class) %>% mutate(x = X - mean(X, na.rm = TRUE),
                                y = Y - mean(Y, na.rm = TRUE),
                                xy = x * y)
  MSD <- heights %>%  dplyr::summarise(MSD = sum( (X - Y) ^2, na.rm = TRUE) / length(X))
  SB <- heights %>%  dplyr::summarise(SB = ( mean(X, na.rm = TRUE) - mean(Y, na.rm = TRUE) )^2 )
  br2 <- heights %>% dplyr::summarise(b = sum(xy, na.rm = TRUE) / sum(x^2, na.rm = TRUE),
                                      r2 = ( sum(xy, na.rm = TRUE)^2 ) / ( sum(x^2, na.rm = TRUE)*sum(y^2, na.rm = TRUE) ))
  br2$NU1 <- (1 - br2$b)^2
  br2$LC1 <- 1 - br2$r2
  NU2 <- heights %>% dplyr::summarise(NU2 = sum(x^2, na.rm = TRUE) / length(x))
  LC2 <- heights %>% dplyr::summarise(LC2 = sum(y^2, na.rm = TRUE) / length(y))
  NU <- cbind(br2, NU2)
  NU$NU <- NU$NU1 * NU$NU2
  LC <- cbind(br2, LC2)
  LC$LC <- LC$LC1 * LC$LC2
  msd <- cbind(MSD, SB)
  msd <- cbind(msd, NU$NU)
  msd <- cbind(msd, LC$LC)
  colnames(msd) <- c('msd', 'sb', 'nu', 'lc')
  msd <- as.data.frame(t(msd))
  colnames(msd) <- 'components'
  msd$var <- rownames(msd)
  msd$landscape <- landscape
  # save for cross landscape comparisons
  write.csv(msd, paste0(evalHeightPath, '/', landscape, '_msd.csv'), row.names = F)

  # # plot msd
  pl2 <- ggplot(msd %>% filter(var != 'msd')) +
  geom_bar(aes(y = components, x = var, fill = var), stat = 'identity') +
  theme_bw()
  pl2


}
