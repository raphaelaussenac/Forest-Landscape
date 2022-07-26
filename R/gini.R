###############################################################
# calculate Gini index on 1ha cells and define stand type
# (even-aged / uneven-aged)
###############################################################

giniClass <- function(tree, df, sce, landscape){

  # load packages
  require(reldist)

  # gini threshold to distinguish even-aged stands from uneven-aged stands
  Gthresh <- 0.45

  # fundtion to define even / uneven stands depending on Gini
  gin <- function(tree, G){
    # if gini < G --> even-aged stand, else --> uneven-aged stand
    gini <- tree %>% group_by(cellID100) %>% mutate(ba = pi * dbh^2 / 4) %>%
                     summarise(gini = gini(x = ba, weights = n),
                               BA = sum((pi * (dbh/200)^2) * n),
                               Dg = sqrt(sum(dbh^2 * n)/sum(n)),
                               meanH = sum(h * n) /sum(n)) %>%
                     mutate(structure = if_else(gini < G, 'even', 'uneven'))
    gini$structure <- as.factor(gini$structure)
    return(gini)
  }


  # define G depending on scenario
  # if scenario is not working for complexity
  if(!('C' %in% sce & landscape == 'bauges')){
    gini <- gin(tree, G = Gthresh)

  # if scenario is working for complexity
  } else if('C' %in% sce & landscape == 'bauges'){

    #  function to count number of even & uneven managed stands
    #  (only managed stands are considered == accessible and not protected)
    countG <- function(df, gin, Gthresh){
      gini <- gin(tree, G = Gthresh)
      gini <- left_join(df, gini, by = 'cellID100') %>% filter(protect == 0, access == 1, !is.na(structure))
      nbEven <- sum(gini$structure == 'even')
      nbUneven <- sum(gini$structure == 'uneven')
      return(c(Gthresh, nbEven, nbUneven))
    }
    nb <- countG(df, gin, Gthresh)

    # if there's more uneven than even stands
    if(nb[2] < nb[3]){
      tab <- data.frame()
      # increase G thrsehold untill reaching a balanced landscape even / uneven
      while(nb[2] < nb[3]){
        Gthresh <- Gthresh + 0.01
        nb <- countG(df, gin, Gthresh)
        print(nb)
        tab <- rbind(tab, nb)
      }
      colnames(tab) <- c('Gthresh', 'nbEven', 'nbUneven')
      tab$diff <- tab$nbEven - tab$nbUneven
      Gthresh <- tab[abs(tab$diff) == min(abs(tab$diff)), 'Gthresh']

    }

    gini <- gin(tree, G = Gthresh)

  }

  return(gini)

}
