heightPred <- function(){

  ###############################################################
  # initialisation
  ###############################################################

  # load packages
  require(nlme)
  require(dplyr)
  require(ggplot2)

  # load models
  load('./data/bauges/hmodels/diam_height_model_18_02_2021.ro')

  # load virtual tree data
  tree <- readRDS(paste0(tempPath, '/trees75.rds'))
  tree <- tree %>% dplyr::select(-wlid)

  # load correspondence between NFI species code and latin Names
  spCor <- read.csv('./data/spCodeCorrespond.csv', sep = ',')

  ###############################################################
  # prepare tree data for predictions
  ###############################################################

  # add salem sp codes
  tree <- merge(tree, spCor[, c('latinName', 'franceCode')], by.x = 'sp', by.y = 'latinName', all.x = TRUE)

  # calculate stand Dg (in 25*25m cells)
  tree <- tree %>% group_by(cellID25) %>%
                   mutate(DgTotFinal = sqrt(sum(dbh^2 * n)/sum(n))) %>%
                   arrange(cellID25) %>%
                   ungroup()
  #
  # calculate relative dbh
  tree$dbh_rel <- tree$dbh / tree$DgTotFinal

  # rename columns
  tree <- tree %>% rename(espar = franceCode) %>% mutate(ipd = 11)

  ###############################################################
  # predict
  ###############################################################

  # predict
  tree$pred <- round(predict(mod_nlme, newdata = tree, level = 0), 2)

  # plot
  ggplot(data = tree[1:10000,], aes(x = dbh, y = pred, col = sp)) +
  geom_point(alpha = 0.5) +
  ylab('h') +
  theme_light()

  # save
  tree <- tree[, c('cellID25', 'sp', 'n', 'dbh', 'pred')] %>%
            rename(h = pred)
  write.csv(tree, paste0(landPath, '/trees.csv'), row.names = F)

}
