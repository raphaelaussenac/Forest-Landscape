heightPred <- function(landscape){

  ###############################################################
  # initialisation
  ###############################################################

  # load packages
  require(nlme)
  require(dtplyr)
  require(dplyr)
  require(ggplot2)

  # load model
  hmodPath <- paste0('./data/', landscape, '/hmodels/')
  roFile <- list.files(path = hmodPath, pattern = '\\.ro$')
  load(paste0(hmodPath, roFile))
  # if landscape == sneznik, also load bauges model to predict
  # height of secondary species
  if(landscape == 'sneznik'){
    saveMod <- mod_nlme
    hmodPath <- paste0('./data/', 'bauges', '/hmodels/')
    roFile <- list.files(path = hmodPath, pattern = '\\.ro$')
    load(paste0(hmodPath, roFile))
    baugesMod <- mod_nlme
    mod_nlme <- saveMod
  }

  # load virtual tree data
  tree <- readRDS(paste0(tempPath, '/trees75.rds'))
  tree <- tree %>% dplyr::select(-wlid)

  # load correspondence between NFI species code and latin Names
  spCor <- read.csv('./data/spCodeCorrespond.csv', sep = ',')
  spCor <- spCor %>% dplyr::select('latinName', 'franceCode') %>% rename(espar = franceCode)

  ###############################################################
  # prepare tree data for predictions
  ###############################################################

  # adapt species names for height models
  if(landscape == 'bauges'){
    tree <- left_join(x = tree, y = spCor, by = c('sp' = 'latinName'))
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
  } else if(landscape == 'sneznik'){
    tree <- left_join(x = tree, y = spCor, by = c('sp' = 'latinName')) %>% rename('esparBauges' = 'espar')
    tree$espar <- NA
    tree[tree$sp == 'Abies alba', 'espar'] <- 'Abal'
    tree[tree$sp == 'Fagus sylvatica', 'espar'] <- 'Fasy'
    tree[tree$sp == 'Picea abies', 'espar'] <- 'Piab'
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

  # function to predict height (parallel ready!)
  predH <- function(landscape, tree, mod_nlme){
    # at sneznik: split into 2 data sets
    # main sp predicted from sneznik h model
    # secondary sp predicted from bauges model
    if(landscape == 'sneznik'){
      # split data
      main <- tree %>% filter(espar != 'OtherSp') %>% dplyr::select(-esparBauges)
      second <- tree %>% filter(espar == 'OtherSp') %>%
                        dplyr::select(-espar) %>%
                        rename('espar' = 'esparBauges')
      # predict
      main$pred <- round(predict(mod_nlme, newdata = main, level = 0), 2)
      second$pred <- round(predict(baugesMod, newdata = second, level = 0), 2)
      # reassemble data
      tree <- bind_rows(main, second)

    } else if(landscape != 'sneznik'){
      # predict
      tree$pred <- round(predict(mod_nlme, newdata = tree, level = 0), 2)
    }
    return(tree)
  }

  tree <- predH(landscape, tree, mod_nlme)

  # save
  tree <- lazy_dt(tree)
  tree <- tree %>% dplyr::select('cellID25', 'sp', 'n', 'dbh', 'pred') %>%
            rename(h = pred) %>% arrange(cellID25) %>% group_by(cellID25) %>%
            arrange(sp, dbh, .by_group = TRUE) %>%
            mutate(across(c(dbh, h), \(x) round(x, 4))) %>%
            ungroup() %>%
            as.data.frame()

  write.csv(tree, paste0(landPath, '/trees.csv'), row.names = F)

  # plot
  ggplot(data = tree[1:30000,], aes(x = dbh, y = h, col = sp)) +
  geom_point(alpha = 0.5) +
  ylab('h') +
  theme_light()

}
