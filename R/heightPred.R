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
tree <- read.csv(paste0(landPath, '/trees75.csv'))

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
tree <- tree %>% rename(espar = franceCode, idp = cellID25)

###############################################################
# predict
###############################################################

# subset tree -> keep only 6 modelled sp and XXX lines
tree <- tree %>% filter(espar %in% c('03', '09', '15S', '17C', '61', '62'))
tree <- tree[1:1000,]

# predict
tree$pred <- round(predict(mod_nlme, newdata = tree, level = 0), 2)






# plot
ggplot(data = tree, aes(x = dbh, y = pred, col = sp)) +
geom_point(alpha = 0.5) +
ylab('h') +
theme_light()




# summary(mod_nlme)


# TODO: what about the other species?
# TODO: how RE are they managed? (level = 0 ok?)
