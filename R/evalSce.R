evalSce <- function(landscape){

  # load package
  require(dplyr)
  require(stringr)
  require(tidyr)

  # load management synthesis tables associated to all scenarios
  managList <- list.files(evalPath, pattern = 'managSyn')

  # function to stack the management tables
  stackT <- function(managT){
    df <- read.csv(paste0(evalPath, '/', managT))
    a <- str_remove_all(managT, pattern = 'managSyn_')
    df$sce <- str_remove_all(a, pattern = '.csv')
    return(df)
  }

  df <- lapply(managList, stackT)
  df <- do.call(rbind, df)

  ###############################################################
  # calculate surfaces
  ###############################################################

  if(landscape == 'bauges'){
    df[is.na(df$manag), 'manag'] <- 'no manag'
    sce <- df %>% dplyr::select(-accessProtec, -owner, -density) %>%
                 mutate(manag = case_when(manag %in% c('coppice', 'final cut') ~ 'even', !(manag %in% c('coppice', 'final cut')) ~ manag)) %>%
                 pivot_longer(!c(manag, sce), names_to = 'sp', values_to = 'count')
  } else if(landscape != 'bauges'){
    sce <- df %>% pivot_longer(!c(manag, sce), names_to = 'sp', values_to = 'count')
  }

  # calculate surfaces

  tot <- sce %>% group_by(sce) %>% summarise(tot = sum(count, na.rm = TRUE))
  unman <- sce %>% group_by(sce, manag) %>% filter(manag == 'no manag') %>% summarise(unman = sum(count, na.rm = TRUE)) %>% dplyr::select(-manag)
  man <- sce %>% group_by(sce) %>% filter(manag != 'no manag') %>% summarise(man = sum(count, na.rm = TRUE))
  even <- sce %>% group_by(sce, manag) %>% filter(manag == 'even') %>% summarise(even = sum(count, na.rm = TRUE)) %>% dplyr::select(-manag)
  uneven <- sce %>% group_by(sce, manag) %>% filter(manag == 'uneven') %>% summarise(uneven = sum(count, na.rm = TRUE)) %>% dplyr::select(-manag)

  # merge tables
  sceSyn <- full_join(tot, unman, by = 'sce') %>%
                full_join(., man, by = 'sce') %>%
                    full_join(., even, by = 'sce') %>%
                        full_join(., uneven, by = 'sce')
  # replace NA with 0
  sceSyn <- mutate_all(sceSyn, ~replace_na(.,0))


  # verifications (should be = 0)
  sceSyn <- sceSyn %>% mutate(totverif = man + unman - tot,
                              manverif = even + uneven - man)

  # calculate proportions
  sceSyn <- sceSyn %>% mutate(propUnman = unman * 100 / tot,
                              propMan = man * 100 / tot,
                              relPropEven = even * 100 / man,
                              relPropUneven = uneven * 100 / man)

  # calculate changes in proportion depending on sce
  baseline <- sceSyn %>% filter(sce == 'B')
  propUnmanB <- baseline$propUnman
  relPropEvenB <- baseline$relPropEven
  relPropUnevenB <- baseline$relPropUneven

  sceSyn <- sceSyn %>% mutate(changeRelPropEven = relPropEven - relPropEvenB,
                              changeRelPropUneven = relPropUneven -relPropUnevenB,
                              changePropUnman = propUnman - propUnmanB)

  # save
  write.csv(sceSyn, paste0(evalPath, '/sceSyn.csv'))
}
