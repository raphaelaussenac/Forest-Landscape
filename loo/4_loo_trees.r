
# load loo results
results <- readRDS(paste0(tempPath, '/loodf_results.rds')) %>% as_tibble()

# load tree data and vegetation type data
tree <- readRDS(paste0(tempPath, '/treeTemp.rds'))
tree$idp <- as.character(tree$idp)

# load packages
library(corrplot)

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
                    rename(BAfield = BAsp.x, BApred = BAsp.y)

rank <- full_join(treeBA, resuBA, by = 'sp') %>%
                    mutate(rank_obs = replace_na(rank_obs, 0),
                           rank_pred = replace_na(rank_pred, 0))


# plot sp abundance
pl1 <- ggplot(data = abun1 %>% filter(sp %in% spList)) +
geom_bar(aes(x = reorder(sp, -BAsp), y = BAsp, fill  = src), stat = 'identity', position = position_dodge()) +
theme_bw() +
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.position = 'bottom', legend.title = element_blank(), axis.title.x = element_blank())

pl2 <- ggplot() +
geom_point(data = abun2, aes(x = BAfield, y = BApred)) +
geom_point(data = abun2 %>% filter(BAfield == 0 | BApred == 0), aes(x = BAfield, y = BApred), col = 'orange') +
theme_bw()

# plot species rank consistency
pl3 <- ggplot() +
geom_point(data = rank, aes(x = rank_pred, y = rank_obs)) +
geom_point(data = rank %>% filter(rank_pred == 0 | rank_obs == 0), aes(x = rank_pred, y = rank_obs), col = 'orange') +
theme_bw()

# pl <- arrangeGrob(pl1, pl2, pl3, ncol = 3, nrow = 1)
# plot(pl)


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

# mod1 <- lm(rank_obs ~ rank_pred, data = rank %>% filter(rank_obs == 1))
# mod2 <- lm(rank_obs ~ rank_pred, data = rank %>% filter(rank_obs == 2))
# mod3 <- lm(rank_obs ~ rank_pred, data = rank %>% filter(rank_obs == 3))
# mod4 <- lm(rank_obs ~ rank_pred, data = rank %>% filter(rank_obs == 4))

pl5 <- ggplot() +
geom_bar(data = rank %>% filter(rank_obs != 999), aes(x = rank_pred, y = rank_obs), stat = 'identity') +
# geom_point(data = rank %>% filter(rank_obs == 1), aes(x = rank_pred, y = rank_obs), col = 'green', alpha = 0.5) +
# geom_point(data = rank %>% filter(rank_obs == 2), aes(x = rank_pred, y = rank_obs), col = 'orange', alpha = 0.5) +
# geom_point(data = rank %>% filter(rank_obs == 3), aes(x = rank_pred, y = rank_obs), col = 'red', alpha = 0.5) +
# geom_point(data = rank %>% filter(rank_obs == 4), aes(x = rank_pred, y = rank_obs), col = 'yellow', alpha = 0.5) +
# coord_fixed() +
facet_wrap(. ~ rank_obs) +
# xlim(0, 90) +
# ylim(0, 90) +
# geom_abline(intercept = 0, slope = 1, color = "black", linetype = 2, size = 1) +
# geom_abline(intercept = mod2$coef[1], slope = mod2$coef[2], color = "red", linetype = 1, size = 1) +
# annotate(geom = 'text', x = 7, y = 3.5, label = paste('R² = ', round(summary(mod2)$r.squared,2)), col = 'red', size = 5) +
#   geom_hex(bins = 8) +
#   scale_fill_continuous(type = "viridis") +
theme_bw()
pl5



pl <- arrangeGrob(pl1, pl2, pl3, pl4, pl5, ncol = 3, nrow = 2)
plot(pl)








# ###############################################################
# # Dominant species (biggest Dg)
# ###############################################################

# # identify species with biggest Dg in field data
# treeDom <- tree %>% group_by(idp, species_name) %>%
#                     summarise(Dg_obs = sqrt(sum(DBH^2 * w)/sum(w))) %>%
#                     mutate(rank_obs = dense_rank(desc(Dg_obs))) %>%
#                     arrange(idp, rank_obs)

# # identify species with biggest Dg in predicted data
# resDom <- results %>% group_by(i, sp) %>%
#                       summarise(Dg_pred = sqrt(sum(dbh^2 * n)/sum(n))) %>%
#                       mutate(rank_pred = dense_rank(desc(Dg_pred))) %>%
#                       arrange(i, rank_pred)


# test <- full_join(treeDom, resDom, join_by('idp' == 'i', 'species_name' == 'sp')) %>%
#                   arrange(idp) %>% na.omit()
#                   mutate(Dg_obs = replace_na(Dg_obs,0),
#                          rank_obs = replace_na(rank_obs,0),
#                          Dg_pred = replace_na(Dg_pred,0),
#                          rank_pred = replace_na(rank_pred,0))
# # test <- na.omit(test)

# M = cor(test[,3:6])
# corrplot(M, method = 'number', type = 'lower') # colorful number



# mod1 <- lm(test$Dg_obs ~ test$Dg_pred)
# summary(mod1)
# plot(test$Dg_obs ~test$Dg_pred)



# mod1 <- lm(test$rank_obs ~ test$rank_pred)
# summary(mod1)
# plot(test$rank_obs ~test$rank_pred)


# ###############################################################
# # most abundant species (highest BA)
# ###############################################################

# # TODO: CorrPlot & hist
# # TODO: territoire et plot!



# # Calculate N, Dg and BA for the whole plote, for deciduous and coniferous and
# # for ech species
# tree <- tree %>% mutate(BAtree =  (pi * (DBH/200)^2) * w) %>%
#                    group_by(idp) %>% mutate(BAtot = sum((pi * (DBH/200)^2) * w),
#                                             Dgtot = sqrt(sum(DBH^2 * w)/sum(w)),
#                                             Ntot = sum(w)) %>%
#                    group_by(idp, spType) %>% mutate(BAdc = sum((pi * (DBH/200)^2) * w),
#                                                     Dgdc = sqrt(sum(DBH^2 * w)/sum(w)),
#                                                     Ndc = sum(w)) %>%
#                    group_by(idp, spType, species_name) %>% mutate(BAsp = sum((pi * (DBH/200)^2) * w),
#                                                                   Dgsp = sqrt(sum(DBH^2 * w)/sum(w))) %>%
#                    group_by(idp, spType) %>% mutate(spPropdc = BAsp/BAdc) %>% ungroup() %>%
#                    mutate(treePropsp = BAtree / BAsp)









# results <- results %>% filter(i == '1007301')
# tree <- tree %>% filter(idp == '1007301')


# ggplot() +
# geom_density(data = results, aes(dbh, weights=n)) +
# geom_density(data = tree, aes(DBH, weights=w), col = 'red')






