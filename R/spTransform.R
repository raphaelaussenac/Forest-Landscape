# transform species code into latin names and create
# deciduous / coniferous categories

spTransform <- function(tree){

  # load correspondence between NFI species code and latin Names
  spCor <- read.csv('./data/spCorrespond.csv', sep = ',')

  # import species latin names
  tree <- merge(tree, spCor[, c('latinName', 'franceCode')], by.x = 'espar', by.y = 'franceCode')
  colnames(tree)[colnames(tree) == 'latinName'] <- 'species_name'

  # correct species name
  tree$species_name <- as.character(tree$species_name)

  # first define deciduous and coniferous species
  deciduousSp <- c('Fagus sylvatica', 'Sorbus aria', 'Quercus robur',
                   'Quercus petraea', 'Acer opalus', 'Tilia platyphyllos',
                   'Carpinus betulus', 'Castanea sativa', 'Populus nigra',
                   'Fraxinus excelsior', 'Corylus avellana', 'Acer campestre',
                   'Betula pendula', 'Populus tremula', 'Prunus avium',
                   'Acer pseudoplatanus', 'Tilia cordata', 'Salix caprea',
                   'Ulmus glabra', 'Sorbus aucuparia', 'Robinia pseudacacia',
                   'Malus sylvestris', 'Acer platanoides', 'Salix alba',
                   'Crataegus monogyna', 'Sorbus mougeoti', 'Laburnum anagyroides',
                   'Quercus pubescens', 'Alnus glutinosa', 'Salix cinerea',
                   'Alnus incana', 'Ilex aquifolium', 'Buxus sempervirens',
                   'Prunus padus', 'Robinia pseudoacacia')
  #
  coniferousSp <- c('Picea abies', 'Abies alba', 'Pinus sylvestris',
                    'Taxus baccata', 'Larix decidua', 'Larix kaempferi',
                    'Picea sitchensis')
  #
  tree[tree$species_name %in% deciduousSp, 'spType'] <- 'D'
  tree[tree$species_name %in% coniferousSp, 'spType'] <- 'C'

  return(tree)

}
