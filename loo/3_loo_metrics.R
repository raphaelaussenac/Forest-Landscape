# load packages
library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyr)

# load loo data
load(paste0('./loo/LOO_values.RDA'))
if(landscape == 'bauges'){
    loodf <- Bauges
    load(paste0('./loo/bauges_tfv.RDA'))
    tfv <- Bauges %>% select(plotMerge, CODE_TFV)
    loodf <- left_join(loodf, tfv, join_by('plot' == 'plotMerge'))
} else if (landscape == 'milicz'){
    loodf <- Milicz %>% mutate(CODE_TFV = 1)
} else if (landscape == 'sneznik'){
    loodf <- Sneznik %>% mutate(CODE_TFV = 1)
    loodf$plot <- as.character(as.numeric(as.factor(loodf$plot)))
}
# loodf <- loodf %>% select(plotMerge, names(Bauges)[str_detect(names(Bauges), 'predicted')], CODE_TFV) %>%
#                    rename(BA = BA_predicted, Dg = Dg_predicted, Dprop = DP_predicted)

# load loo results
results <- readRDS(paste0(tempPath, '/loodf_results.rds')) %>% as_tibble()

# load list of deciduous and coniferous sp
deciduousSp <- readRDS('./data/deciduousSp.rds')
coniferousSp <- readRDS('./data/coniferousSp.rds')

###############################################################
# convert DP_field DP_predicted into %
###############################################################

loodf$DP_predicted <- loodf$DP_predicted * 100
loodf$DP_field <- loodf$DP_field * 100


###############################################################
# calculate BA, Dg & Dprop on results
###############################################################

# define spType
results[results$sp %in% deciduousSp, 'spType'] <- 'D'
results[results$sp %in% coniferousSp, 'spType'] <- 'C'
results$spType <- as.factor(results$spType)

# Calculate N, Dg and BA for the whole plote, for deciduous and coniferous and
# for ech species
df <- results %>% group_by(i) %>% summarise(BAtot = 16 * sum((pi * (dbh/200)^2) * n),
                                            Dgtot = sqrt(sum(dbh^2 * n)/sum(n)))

dc <- results %>% group_by(i, spType) %>% summarise(BAtot = 16 * sum((pi * (dbh/200)^2) * n))
dc <- dc %>% pivot_wider(names_from = spType, values_from = BAtot) %>%
             mutate(D = replace_na(D,0),
                    C = replace_na(C,0),
                    Dprop = (D / (D+C)) * 100) %>% select(i, Dprop)

# merge results with loodf
df <- left_join(loodf, df, join_by('plot' == 'i'))
df <- df %>% filter(!is.na(CODE_TFV))
df <- left_join(df, dc, join_by('plot' == 'i'))



###############################################################
# plot
###############################################################

dfBackup <- df

# total basal area
df <- dfBackup
df <- dfBackup %>% filter(!is.na(BA_field), !is.na(BAtot))
mod1 <- lm(df$BA_field ~ df$BAtot)
rmse1 <- sqrt(mean(mod1$residuals^2))
rmse11 <- sqrt( sum( (df$BAtot - df$BA_field)^2 ) / nrow(df) )
rmsre1 <- sqrt( sum( (100 - ( (df$BAtot * 100) / df$BA_field ))^2 ) / nrow(df) )

pl1 <- ggplot(data = df) +
geom_point(aes(x = BAtot, y = BA_field), pch = 16) +
coord_fixed() +
xlim(0, 90) +
ylim(0, 90) +
labs(y = expression(BA[field])) +
labs(x = expression(BA[pred])) +
geom_abline(intercept = 0, slope = 1, color = "black", linetype = 2, size = 1) +
# geom_abline(intercept = mod1$coef[1], slope = mod1$coef[2], color = "red", linetype = 1, size = 1) +
# annotate(geom = 'text', x = 75, y = 40, label = paste('R² = ', round(summary(mod1)$r.squared,2)), col = 'red', size = 5) +
# annotate(geom = 'text', x = 75, y = 35, label = paste('RMSE = ', round(rmse1,2)), col = 'red', size = 5) +
annotate(geom = 'text', x = 75, y = 30, label = paste('RMSE = ', round(rmse11,2)), col = 'red', size = 5) +
annotate(geom = 'text', x = 75, y = 20, label = paste('RMSRE = ', round(rmsre1,2)), col = 'red', size = 5) +
theme_classic() 


# quadratic diameter
df <- dfBackup
df <- dfBackup %>% filter(!is.na(Dg_field), !is.na(Dgtot))
mod2 <- lm(df$Dg_field ~ df$Dgtot)
rmse2 <- sqrt(mean(mod2$residuals^2))
rmse22 <- sqrt( sum( (df$Dgtot - df$Dg_field)^2 ) / nrow(df) )
rmsre2 <- sqrt( sum( (100 - ( (df$Dgtot * 100) / df$Dg_field ))^2 ) / nrow(df) )

pl2 <- ggplot(data = df) +
geom_point(aes(x = Dgtot, y = Dg_field), pch = 16) +
coord_fixed() +
xlim(5, 100) +
ylim(5, 100) +
labs(y = expression(Dg[field])) +
labs(x = expression(Dg[pred])) +
geom_abline(intercept = 0, slope = 1, color = "black", linetype = 2, size = 1) +
# geom_abline(intercept = mod2$coef[1], slope = mod2$coef[2], color = "red", linetype = 1, size = 1) +
# annotate(geom = 'text', x = 70, y = 40, label = paste('R² = ', round(summary(mod2)$r.squared,2)), col = 'red', size = 5) +
# annotate(geom = 'text', x = 70, y = 35, label = paste('RMSE = ', round(rmse2,2)), col = 'red', size = 5) +
annotate(geom = 'text', x = 70, y = 30, label = paste('RMSE = ', round(rmse22,2)), col = 'red', size = 5) +
annotate(geom = 'text', x = 70, y = 20, label = paste('RMSRE = ', round(rmsre2,2)), col = 'red', size = 5) +
theme_classic() 


# deciduous proportion
df <- dfBackup
df <- dfBackup %>% filter(!is.na(DP_field), !is.na(Dprop))
mod3 <- lm(df$DP_field ~ df$Dprop)
rmse3 <- sqrt(mean(mod3$residuals^2))
rmse33 <- sqrt( sum( (df$Dprop - df$DP_field)^2 ) / nrow(df) )
rmsre3 <- sqrt( sum( (100 - ( (df$Dprop * 100) / df$DP_field ))^2 ) / nrow(df) )


pl3 <- ggplot(data = df) +
geom_point(aes(x = Dprop, y = DP_field), pch = 16) +
coord_fixed() +
# xlim(10, 80) +
ylim(0, 100) +
labs(y = expression(Dprop[field])) +
labs(x = expression(Dprop[pred])) +
geom_abline(intercept = 0, slope = 1, color = "black", linetype = 2, size = 1) +
# geom_abline(intercept = mod3$coef[1], slope = mod3$coef[2], color = "red", linetype = 1, size = 1) +
# annotate(geom = 'text', x = 80, y = 55, label = paste('R² = ', round(summary(mod3)$r.squared,2)), col = 'red', size = 5) +
# annotate(geom = 'text', x = 80, y = 50, label = paste('RMSE = ', round(rmse3,2)), col = 'red', size = 5) +
annotate(geom = 'text', x = 80, y = 45, label = paste('RMSE = ', round(rmse33,2)), col = 'red', size = 5) +
annotate(geom = 'text', x = 80, y = 35, label = paste('RMSRE = ', round(rmsre3,2)), col = 'red', size = 5) +
theme_classic() 

pl <- arrangeGrob(pl1, pl2, pl3, ncol = 3, nrow = 1)
plot(pl)
saveRDS(pl, paste0('./loo/', landscape, '_loo.RDS'))

# ggsave(file = 'loo.pdf', plot = pl, width = 6, height = 16)


###############################################################
# 
###############################################################

# par(mfrow=c(1, 3))

# plot(df$BA_predicted ~ df$BAtot, ylab = 'BA lidar pred', xlab = 'BA virtual tree pred')

# plot(df$Dg_predicted ~ df$Dgtot, ylab = 'Dg lidar pred', xlab = 'Dg virtual tree')

# plot(df$DP_predicted ~ df$Dprop, ylab = 'Dprop lidar pred', xlab = 'Dprop virtual tree')

