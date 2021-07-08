# transform species code into latin names and create
# deciduous / coniferous categories

spTransform <- function(tree){

  # load correspondence between NFI species code and latin Names
  spCor <- read.csv('./data/spCodeCorrespond.csv', sep = ',')

  # import species latin names
  tree <- merge(tree, spCor[, c('latinName', 'franceCode')], by.x = 'espar', by.y = 'franceCode', all.x = TRUE)
  colnames(tree)[colnames(tree) == 'latinName'] <- 'species_name'

  # correct species name
  tree$species_name <- as.character(tree$species_name)

  # load list of deciduous and coniferous sp
  deciduousSp <- readRDS('./data/deciduousSp.rds')
  coniferousSp <- readRDS('./data/coniferousSp.rds')
  
  tree[tree$species_name %in% deciduousSp, 'spType'] <- 'D'
  tree[tree$species_name %in% coniferousSp, 'spType'] <- 'C'

  return(tree)

}
