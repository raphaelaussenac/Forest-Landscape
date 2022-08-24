evalHeight <- function(landscape){

  ###############################################################
  # initialisation
  ###############################################################

  # load packages
  require(raster)
  require(dplyr)
  require(ggplot2)
  require(ggExtra)
  require(stringr)

  # load LiDAR heights
  lidH <- list.files(paste0('./data/', landscape, '/LiDARh'), pattern = '.rda')
  load(paste0('./data/', landscape, '/LiDARh/', lidH))

  # load tree data
  tree <- read.csv(paste0(landPath, '/trees.csv'))

  # load cellID25 raster
  cellID25 <- raster(paste0(landPath, '/cellID25.asc'))
  names(cellID25) <- 'cellID25'

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
    compH <- tree %>% group_by(cellID25) %>%
                      arrange(-cellID25, -h) %>%
                      mutate(cumN = cumsum(n)) %>%
                      mutate(NN = ifelse(cumN<=nbTrees, 0, cumN-nbTrees)) %>%
                      mutate(n2 = abs(n - NN)) %>%
                      mutate(cumN2 = cumsum(n2)) %>%
                      filter(cumN2 <= nbTrees) %>%
                      select(-n, cumN, -NN, -cumN2) %>%
                      summarise(h = weighted.mean(h, n2), hlid = unique(get(paste0('Hm', nbTrees))))
    return(compH)
  }
  heights <- compareH(6, tree)
  mod <- lm(hlid~h, data = heights)
  summary(mod)

  # create classes for h --> boxplot
  heights <- heights %>% mutate(bin=cut_width(h, width = 1, center = 0)) %>%
                         group_by(bin) %>% mutate(class = round(mean(h))) %>% ungroup()

  # plot
  pl1 <- ggplot(heights, aes(x = h, y = hlid)) +
  geom_point(alpha = 0, size = 0.5, col = 'grey', pch = 16) +
  xlim(5, 45) +
  ylim(5, 45) +
  labs(y = expression(Hdom[L])) +
  labs(x = expression(Hdom[T])) +
  coord_fixed() +
  annotate(geom = 'text', x = 15, y = 40, label = paste('R² = ', round(summary(mod)$r.squared,2)), col = 'red', size = 5) +
  ggtitle(str_to_title(landscape)) +
  geom_boxplot(aes(x = class, y = hlid, group = class), alpha = 0.8, outlier.shape = NA) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = 2, size = 1) +
  geom_abline(intercept = mod$coef[1], slope = mod$coef[2], color = "red", linetype = 1, size = 1) +
  # theme_bw() +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))
  pl1 <- ggMarginal(pl1, type = 'histogram', margin = 'x', size = 7, xparams = list(fill = "white", bins = 100))
  # pl1
  ggsave(file = paste0(evalHeightPath, '/H', landscape, '.pdf'), plot = pl1, width = 7, height = 7)

  # save ggplot object
  saveRDS(pl1 , paste0('./evalHeight/heights', '_', landscape, '.rds'))

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
