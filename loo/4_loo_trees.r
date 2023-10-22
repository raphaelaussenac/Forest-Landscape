# load packages
library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyr)
library(openxlsx)
library(Hmisc)
library(nlme)
library(stringr)

# function to replace commas by dots and transform character vectors
  # into numeric vectors
  comTodot <- function(charVect){
    charVect <- as.numeric(str_replace_all(charVect, ',', '\\.'))
    return(charVect)
  }



if(landscape == 'bauges'){
    # load loo results
    results <- readRDS(paste0(tempPath, '/loodf_results.rds')) %>%
                        filter(dbh>17.5) %>%
                        as_tibble()
    # load ONF sp codes
    spCodes <- read.xlsx('./loo/dico_attribut_latin.xlsx', sheet = "liste_essence") %>%
                     select(Cod_ess, Libelle_ess)
    latin <- read.xlsx('./loo/dico_attribut_latin.xlsx', sheet = "latin")

    # load tree data and vegetation type data
    load('./loo/trees_Bauges_328_plots_with_duplicates.RDA')
    tree <- trees_Bauges %>%
            select(-Fam_ess2, -Fam_ess3, -plotID) %>%
            rename(idp = plotMerge, DBH = Diam, spType  = Fam_ess1, species_name = Cod_ess) %>%
            mutate(w = 14.14711) # plot radius = 15. pi*15^2 = 706.8583. weight = 10000/(pi*15^2) = 14.14711
#     tree <- left_join(tree, tfv, join_by(idp == plotMerge))
#     tree <- groupTfv(tree) %>%
#             relocate(idp, w, CODE_TFV, species_name, spType, DBH)
    tree <- left_join(tree, spCodes, join_by(species_name == Cod_ess))
    tree <- left_join(tree, latin, join_by(Libelle_ess == fr)) %>%
            select(-species_name, - Libelle_ess) %>%
            rename(species_name = latin)
} else if (landscape == 'milicz'){
    # load loo results
    results <- readRDS(paste0(tempPath, '/loodf_results.rds')) %>% as_tibble()
    # load tree data and vegetation type data
    tree <- readRDS(paste0(tempPath, '/treeTemp.rds'))
    tree$idp <- as.character(tree$idp)
} else if (landscape == 'sneznik'){
    # load loo results
    results <- readRDS(paste0(tempPath, '/loodf_results.rds')) %>% as_tibble()
    # load tree data and vegetation type data
    tree <- readRDS(paste0(tempPath, '/treeTemp.rds'))
    tree$idp <- as.numeric(as.factor(tree$idp))
    tree$idp <- as.character(tree$idp)
}


###############################################################
# dbh quantile validation
###############################################################

# get dbh quantile 0.95
quant <- 0.95
qtree <- tree %>% group_by(idp) %>% summarise(obsq95 = wtd.quantile(DBH, weights = w, probs = quant))
qresults <- results %>% group_by(i) %>% summarise(predq95 = wtd.quantile(dbh, weights = n, probs = quant))
q <- inner_join(qtree, qresults, join_by(idp == i))

# plot
mod1 <- lm(obsq95 ~ predq95, data = q)
rmse <- sqrt( sum( (q$predq95 - q$obsq95)^2 ) / nrow(q) )
rmsre <- sqrt( sum( (100 - ( (q$predq95 * 100) / q$obsq95 ))^2 ) / nrow(df) )


dbhQ <- ggplot(data = q) +
geom_point(aes(x = predq95, y = obsq95), pch = 16) +
coord_fixed() +
# xlim(0, 110) +
# ylim(0, 110) +
labs(y = expression(DBH(Q[0.95])[field])) +
labs(x = expression(DBH(Q[0.95])[pred])) +
geom_abline(intercept = 0, slope = 1, color = "black", linetype = 2, size = 1) +
# geom_abline(intercept = mod1$coef[1], slope = mod1$coef[2], color = "red", linetype = 1, size = 1) +
# annotate(geom = 'text', x = 75, y = 40, label = paste('R² = ', round(summary(mod1)$r.squared,2)), col = 'red', size = 5) +
# annotate(geom = 'text', x = 75, y = 35, label = paste('RMSE = ', round(rmse1,2)), col = 'red', size = 5) +
annotate(geom = 'text', x = 75, y = 30, label = paste('RMSE = ', round(rmse,2)), col = 'red', size = 5) +
annotate(geom = 'text', x = 75, y = 20, label = paste('RMSRE = ', round(rmsre,2)), col = 'red', size = 5) +
theme_classic()

saveRDS(dbhQ, paste0('./loo/', landscape, '_dbhQ.RDS'))

# x <- c(6, 6, 6, 2, 2, 2, 3, 4, 5, 5, 5, 5, 1, 1, 1, 1, 1, 1, 1)
# w <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
# quantile(x, probs=c(0, .25, .5, .75, 1))
# wtd.quantile(x, weights=w, probs=c(0, .25, .5, .75, 1))#,

# x <- c(6, 2, 3, 4, 5, 1)
# w <- c(3, 3, 1, 1, 4, 7)
# wtd.quantile(x, weights=w, probs=c(0, .25, .5, .75, 1))#,

###############################################################
# h quantile validation
###############################################################


if (landscape != 'bauges'){
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

       # load correspondence between NFI species code and latin Names
       spCor <- read.csv('./data/spCodeCorrespond.csv', sep = ',')
       spCor <- spCor %>% dplyr::select('latinName', 'franceCode') %>% rename(espar = franceCode)

       # adapt species names for height models
       if(landscape == 'milicz'){
       results$espar <- NA
       results[results$sp == 'Pinus sylvestris', 'espar'] <- 'Pisy'
       results[results$sp == 'Fagus sylvatica', 'espar'] <- 'Fasy'
       results[results$sp == 'Picea abies', 'espar'] <- 'Piab'
       results[results$sp == 'Quercus robur', 'espar'] <- 'Quun'
       results[results$sp == 'Betula pendula', 'espar'] <- 'Bepe'
       results[results$sp == 'Alnus glutinosa', 'espar'] <- 'Algl'
       results[results$sp == 'Carpinus betulus', 'espar'] <- 'Cabe'
       results[results$sp == 'Larix decidua', 'espar'] <- 'Lade'
       results[results$sp == 'Tilia cordata', 'espar'] <- 'Tico'
       results[results$sp == 'Quercus rubra', 'espar'] <- 'Quru'
       results[results$sp == 'Acer pseudoplatanus', 'espar'] <- 'Acps'
       results[results$sp == 'Prunus serotina', 'espar'] <- 'Prse'
       results[is.na(results$espar), 'espar'] <- 'OtherSp'
       } else if(landscape == 'sneznik'){
       results <- left_join(x = results, y = spCor, by = c('sp' = 'latinName')) %>% rename('esparBauges' = 'espar')
       results$espar <- NA
       results[results$sp == 'Abies alba', 'espar'] <- 'Abal'
       results[results$sp == 'Fagus sylvatica', 'espar'] <- 'Fasy'
       results[results$sp == 'Picea abies', 'espar'] <- 'Piab'
       results[is.na(results$espar), 'espar'] <- 'OtherSp'
       }

       # calculate stand Dg (in 25*25m cells)
       h <- results %>% group_by(i) %>%
                     mutate(DgTot = sqrt(sum(dbh^2 * n)/sum(n))) %>%
                     arrange(i) %>%
                     ungroup()

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

       h <- predH(landscape, h, mod_nlme)


       # get heights of trees on the field
       if (landscape == 'milicz'){
              fieldH <- read.csv('./data/milicz/inventory/Milicz_trees.csv', sep = ';')
              fieldH$DBH <- comTodot(fieldH$DBH)
              fieldH$height <- comTodot(fieldH$height)
              fieldH$species <- as.factor(fieldH$species)
              fieldH$species <- recode_factor(fieldH$species, 'Malus silvestris' = 'Malus sylvestris',
                                                               'Quercus undefined' = 'Quercus robur',
                                                               'Rhamnus frangula' = 'Frangula alnus',
                                                               'Ulmus' = 'Ulmus minor')

              # add missing columns remove useless columns
              fieldH <- fieldH %>% dplyr::select(plot_id, species, height) %>%
                            rename(i = plot_id, sp = species, Hfield = height) %>%
                            mutate(i = as.factor(i), n = NA) %>%
                            relocate(i, n, sp)
              #
              # define fieldH weights
              # plot radius = 12.62 m --> pi*12.62^2 = 500.3439 m2
              fieldH$n <- 10000 / 500.3439

       } else if (landscape == 'sneznik'){
              # load tree inventory data
              fieldH <- read.xlsx('./data/sneznik/inventory/Sneznik_tree_data.xlsx', sheet = 1)

              # replace plot id with letters by numeric ids
              fieldH$PLOTID <- as.numeric(as.factor(fieldH$PLOTID))

              # rename undefined species
              fieldH$SPECIES <- as.factor(fieldH$SPECIES)
              fieldH$SPECIES <- recode_factor(fieldH$SPECIES, 'Tilia sp.' = 'Tilia platyphyllos',
                                                        'Salix sp.' = 'Salix caprea',
                                                        'Pinus sp.' = 'Pinus sylvestris')

              # add missing columns remove useless columns
              fieldH <- fieldH %>% dplyr::select(PLOTID, H, DBH, SPECIES) %>%
                            rename(i = PLOTID, sp = SPECIES, Hfield = H) %>%
                            mutate(i = as.factor(i), n = NA) %>%
                            relocate(i, n, sp, DBH)
              #
              # define fieldH weights
              fieldH <- fieldH %>% mutate(n = case_when(DBH <= 29.9 ~ 10000 / (pi * 7.98^2), DBH > 29.9 ~ 10000 / (pi * 12.62^2))) %>%
                                   filter(Hfield > 0)
              # number of H measurments / plot
              test <- fieldH %>% group_by(i) %>% summarise(c = n())
              table(test$c)

       }


# get h quantile 0.95
quant <- 0.95
fieldH <- fieldH %>% group_by(i) %>% summarise(obsq95 = wtd.quantile(Hfield, weights = n, probs = quant))
h <- h %>% group_by(i) %>% summarise(predq95 = wtd.quantile(pred, weights = n, probs = quant))
h <- inner_join(h, fieldH, by = 'i')

# plot
modh <- lm(obsq95 ~ predq95, data = h)
rmse <- sqrt( sum( (h$predq95 - h$obsq95)^2 ) / nrow(h) )
rmsre <- sqrt( sum( (100 - ( (h$predq95 * 100) / h$obsq95 ))^2 ) / nrow(df) )

HQ <- ggplot(data = h) +
geom_point(aes(x = predq95, y = obsq95), pch = 16) +
coord_fixed() +
# xlim(0, 50) +
# ylim(0, 50) +
labs(y = expression(H(Q[0.95])[field])) +
labs(x = expression(H(Q[0.95])[pred])) +
geom_abline(intercept = 0, slope = 1, color = "black", linetype = 2, size = 1) +
# geom_abline(intercept = modh$coef[1], slope = modh$coef[2], color = "red", linetype = 1, size = 1) +
# annotate(geom = 'text', x = 35, y = 20, label = paste('R² = ', round(summary(modh)$r.squared,2)), col = 'red', size = 5) +
# annotate(geom = 'text', x = 75, y = 35, label = paste('RMSE = ', round(rmse1,2)), col = 'red', size = 5) +
annotate(geom = 'text', x = 20, y = 25, label = paste('RMSE = ', round(rmse,2)), col = 'red', size = 5) +
annotate(geom = 'text', x = 20, y = 15, label = paste('RMSRE = ', round(rmsre,2)), col = 'red', size = 5) +
theme_classic()


saveRDS(HQ, paste0('./loo/', landscape, '_', quant, '_hQ.RDS'))

}


###############################################################
# Sp abundance at the landscape level
###############################################################

# Species BA over landscape
treeBA <- tree %>% group_by(species_name) %>%
                   reframe(BAsp = sum((pi * (DBH/200)^2) * w)) %>%
                   arrange(-BAsp) %>%
                   mutate(src = 'field') %>%
                   mutate(rank_obs = dense_rank(desc(BAsp))) %>%
                   rename(sp = species_name) %>% as.data.frame()
spList <- treeBA[1:10, 'sp']

resuBA <- results %>% group_by(sp) %>%
                   reframe(BAsp = 16 * sum((pi * (dbh/200)^2) * n)) %>%
                   arrange(-BAsp) %>%
                   mutate(src = 'pred') %>%
                   mutate(rank_pred = dense_rank(desc(BAsp))) 

abun1 <- bind_rows(treeBA, resuBA) %>% arrange(-BAsp)
abun2 <- full_join(treeBA, resuBA, by = 'sp') %>%
                    mutate(BAsp.x = replace_na(BAsp.x, 0),
                           BAsp.y = replace_na(BAsp.y, 0)) %>%
                    rename(BAfield = BAsp.x, BApred = BAsp.y) %>%
                    arrange(-BAfield) %>%
                    mutate(BAtot = sum(BAfield),
                           BAprop = BAfield * 100 / BAtot,
                           BApropcum = cumsum(BAprop))
rarespList <- abun2 %>% filter(BAprop < 1) %>% select(sp) %>% as.data.frame()
rarespList <- rarespList$sp

rank <- full_join(treeBA, resuBA, by = 'sp') %>%
                    mutate(rank_obs = replace_na(rank_obs, 0),
                           rank_pred = replace_na(rank_pred, 0))


# plot sp abundance
diff <- abun2 %>% mutate(reldiff = round((BApred * 100 /BAfield)-100, 2),
                         diff = round(BApred - BAfield, 2))
abun1 <- left_join(abun1, diff %>% select(sp, diff, reldiff), join_by(sp)) %>%
         filter(sp %in% spList) %>%
         mutate(src = case_when(src == 'field' ~ 'observation',
                                src == 'pred' ~ 'prediction'))
abun12 <- abun1 %>% filter(src == 'prediction')
head(abun12)

pl1 <- ggplot() +
geom_bar(data = abun1, aes(x = reorder(sp, -BAsp), y = BAsp, fill  = src), stat = 'identity', position = position_dodge()) +
theme_classic() +
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.position = 'none', legend.title = element_blank(), axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
# geom_text(data = abun12, aes(x = sp, y = BAsp, label = reldiff), vjust = -0.3, size = 3.5, col = 'black') +
labs(y = 'BA') +
ggtitle(str_to_title(landscape))
pl1
saveRDS(pl1, paste0('./loo/', landscape, '_sp_compo.RDS'))

pl2 <- ggplot() +
geom_point(data = abun2, aes(x = BAfield, y = BApred)) +
geom_point(data = abun2 %>% filter(BAfield == 0 | BApred == 0), aes(x = BAfield, y = BApred), col = 'orange') +
theme_bw()

# plot species rank consistency
pl3 <- ggplot() +
geom_point(data = rank %>% filter(!(sp %in% rarespList)), aes(x = rank_pred, y = rank_obs), size = 2) +
geom_point(data = rank %>% filter(sp %in% rarespList), aes(x = rank_pred, y = rank_obs), col = 'blue', size = 1) +
geom_point(data = rank %>% filter(rank_pred == 0 | rank_obs == 0), aes(x = rank_pred, y = rank_obs), col = 'red', shape = 4, size = 2) +
theme_bw()
pl3

###############################################################
# Sp abundance at the plot level
###############################################################

# Species BA over landscape
treeBA <- tree %>% group_by(species_name, idp) %>%
                   summarise(BAsp_obs = sum((pi * (DBH/200)^2) * w)) %>%
                   mutate(src = 'tree') %>%
                   group_by(idp) %>%
                   mutate(rank_obs = dense_rank(desc(BAsp_obs))) %>%
                   rename(sp = species_name, i = idp) %>% as.data.frame() %>%
                   arrange(i, rank_obs)
# spList <- treeBA[1:10, 'sp']

resuBA <- results %>% group_by(sp, i) %>%
                   summarise(BAsp_pred = 16 * sum((pi * (dbh/200)^2) * n)) %>%
                   group_by(i) %>%
                   mutate(rank_pred = dense_rank(desc(BAsp_pred))) %>%
                   mutate(src = 'resu')%>%
                   arrange(i, rank_pred)

abun2 <- full_join(treeBA, resuBA, join_by(sp == sp, i == i)) %>%
                    mutate(BAsp_obs = replace_na(BAsp_obs, 0),
                           BAsp_pred = replace_na(BAsp_pred, 0))
rank <- full_join(treeBA, resuBA, join_by(sp == sp, i == i)) %>%
                    mutate(rank_obs = replace_na(rank_obs, 999),
                           rank_pred = replace_na(rank_pred, 999)) %>%
                    mutate(rank_pred = as.character(rank_pred),
                           rank_pred = case_when(rank_pred == '999' ~ 'NP',
                                                 TRUE ~ rank_pred))


# plot sp abundance
mod1 <- lm(BAsp_obs ~ BAsp_pred, data = abun2 %>% filter(rank_obs == 1))
mod2 <- lm(BAsp_obs ~ BAsp_pred, data = abun2 %>% filter(rank_obs == 2))
mod3 <- lm(BAsp_obs ~ BAsp_pred, data = abun2 %>% filter(rank_obs == 3))
mod4 <- lm(BAsp_obs ~ BAsp_pred, data = abun2 %>% filter(rank_obs == 4))

pl4 <- ggplot() +
geom_point(data = abun2, aes(x = BAsp_pred, y = BAsp_obs), col = 'green', alpha = 0.5) +
geom_point(data = abun2 %>% filter(rank_obs == 1), aes(x = BAsp_pred, y = BAsp_obs), col = 'green', alpha = 0.5) +
geom_point(data = abun2 %>% filter(rank_obs == 2), aes(x = BAsp_pred, y = BAsp_obs), col = 'orange', alpha = 0.5) +
geom_point(data = abun2 %>% filter(rank_obs == 3), aes(x = BAsp_pred, y = BAsp_obs), col = 'red', alpha = 0.5) +
geom_point(data = abun2 %>% filter(rank_obs == 4), aes(x = BAsp_pred, y = BAsp_obs), col = 'yellow', alpha = 0.5) +
coord_fixed() +
# xlim(0, 90) +
# ylim(0, 90) +
geom_abline(intercept = 0, slope = 1, color = "black", linetype = 2, size = 1) +
geom_abline(intercept = mod1$coef[1], slope = mod1$coef[2], color = "green4", linetype = 1, size = 1) +
geom_abline(intercept = mod2$coef[1], slope = mod2$coef[2], color = "orange4", linetype = 1, size = 1) +
geom_abline(intercept = mod3$coef[1], slope = mod3$coef[2], color = "red4", linetype = 1, size = 1) +
geom_abline(intercept = mod4$coef[1], slope = mod4$coef[2], color = "yellow4", linetype = 1, size = 1) +
# annotate(geom = 'text', x = 45, y = 25, label = paste('R² = ', round(summary(mod1)$r.squared,2)), col = 'green4', size = 5) +
theme_bw()
pl4





# plot species rank consistency
rank2 <- rank %>% filter(rank_obs != 999) %>%
              #    select(-BAsp_obs, -src.x, -BAsp_pred, -src.y) %>%
                 group_by(rank_obs) %>%
              #    arrange(rank_obs) %>% 
                 count(rank_pred) %>%
                 mutate(n_tot = sum(n),
                        n_prop = round(n * 100 / n_tot, 2))
rank2 <- rank2 %>% filter(rank_obs < 5)
pl5 <- ggplot() +
geom_bar(data = rank2, aes(x = rank_pred, y = n), stat = 'identity') +
# geom_bar(data = rank %>% filter(rank_obs != 999), aes(x = rank_pred, y = rank_obs), stat = 'identity') +
# ylim(0, 500) +
facet_wrap(. ~ rank_obs) +
# annotate(geom = 'text', x = 75, y = 35, label = paste('RMSE = ', round(rmse1,2)), col = 'red', size = 5) +
# geom_text(aes(x = CODE_TFV, y = N * 20, label = N), vjust = -0.3, size = 3.5, col = 'orange') +
# geom_text(data = surfArea[3,], aes(x = CODE_TFV, y = surface, label =  paste(round(prop, 2), '% of forest cover')), vjust = -1, hjust = 0, size = 3.5, col = 'black') +
geom_text(data = rank2, aes(x = rank_pred, y = n, label = n_prop), vjust = -0.3, size = 3.5, col = 'orange') +
theme_bw()
pl5



# pl <- arrangeGrob(pl1, pl2, pl3, pl4, pl5, ncol = 3, nrow = 2)
# plot(pl)

pl <- arrangeGrob(pl1, pl3, pl5, ncol = 3, nrow = 1)
plot(pl)
saveRDS(pl, paste0('./loo/', landscape, '_compo.RDS'))
# ggsave(file = './loo/compo.pdf', plot = pl, width = 20, height = 20)

