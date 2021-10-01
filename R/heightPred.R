heightPred <- function(landscape){

  ###############################################################
  # initialisation
  ###############################################################

  # load packages
  require(nlme)
  require(dplyr)
  require(ggplot2)

  # load model
  hmodPath <- paste0('./data/', landscape, '/hmodels/')
  roFile <- list.files(path = hmodPath, pattern = '\\.ro$')
  load(paste0(hmodPath, roFile))

  # load virtual tree data
  tree <- readRDS(paste0(tempPath, '/trees75.rds'))
  tree <- tree %>% dplyr::select(-wlid)

  # load correspondence between NFI species code and latin Names
  spCor <- read.csv('./data/spCodeCorrespond.csv', sep = ',')

  ###############################################################
  # prepare tree data for predictions
  ###############################################################

  # adapt species names for height models
  if(landscape == 'bauges'){
    tree <- merge(tree, spCor[, c('latinName', 'franceCode')], by.x = 'sp', by.y = 'latinName', all.x = TRUE)
    tree <- tree %>% rename(espar = franceCode)
  } else if(landscape == 'milicz'){
    tree$espar <- NA
    tree[tree$sp == 'Pinus sylvestris', 'espar'] <- 'Pisy'
    tree[tree$sp == 'Fagus sylvatica', 'espar'] <- 'Fasy'
    tree[tree$sp == 'Picea abies', 'espar'] <- 'Piab'
    tree[tree$sp == 'Quercus robur', 'espar'] <- 'Quun'
    tree[tree$sp == 'Betula pendula', 'espar'] <- 'Bepe'
    tree[tree$sp == 'Alnus glutinosa', 'espar'] <- 'Algl'
    tree[tree$sp == 'Carpinus betulus', 'espar'] <- 'Cabe'
    tree[tree$sp == 'Larix decidua', 'espar'] <- 'Lade'
    tree[tree$sp == 'Tilia cordata', 'espar'] <- 'Tico'
    tree[tree$sp == 'Quercus rubra', 'espar'] <- 'Quru'
    tree[tree$sp == 'Acer pseudoplatanus', 'espar'] <- 'Acps'
    tree[tree$sp == 'Prunus serotina', 'espar'] <- 'Prse'
    tree[is.na(tree$espar), 'espar'] <- 'OtherSp'
  }

  # calculate stand Dg (in 25*25m cells)
  tree <- tree %>% group_by(cellID25) %>%
                   mutate(DgTot = sqrt(sum(dbh^2 * n)/sum(n))) %>%
                   arrange(cellID25) %>%
                   ungroup()


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
