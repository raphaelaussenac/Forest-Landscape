# variable selection algo
#
# 1 - all possible combination of variables are compared with the AIC
# 2 - the model with the largest number of significant terms is chosen
#     among the models with AIC < AICmin + 2
# 3 - linear effects are manually added (if missing) when quadratic term are
#     present. Terms becoming non-signifant because of this operation are removed
#     until all remaining terms are significant
# --------------------------------------------
# 4 - The beech model showed hardly interpretable effects of expoNS and expoEW
#     (+ their quadratic effects). This may have been due to the correlation of
#     their coef in the model (with each other and with slope and slope2). Some
#     of these correlation where > 0.5 (a threshold not found elsewhere in the
#     other models). Removing one of the expo quadratic effects make the other
#     expo terms non-significant, which led us to remove all 4 expo terms.


###############################################################
# initialisation
###############################################################

# clean up environment
rm(list = ls())

# load packages
library(rgdal)
library(raster)
library(sf)
library(dplyr)
library(glmulti)
library(ggplot2)

# set work directory
setwd('C:/Users/raphael.aussenac/Documents/GitHub/LandscapeInit')

# load species-specific SI and env variables at NFI plots
salemSI <- read.csv('./data/salemSI/bdBauges_for_SI_calibration_2021_03_11.txt', sep = '\t')

###############################################################
# manage data format, variable class, unit, model selection...
###############################################################

# remove useless variables
salemSI <- salemSI %>% select(-idp, -xl93, -yl93, -aspect)

# define variable class
salemSI$GRECO <- as.factor(salemSI$GRECO)

# rename variable
salemSI <- salemSI %>% rename(pH = pH_decor)

# add quadratic effects
salemSI <- salemSI %>% mutate(elev2 = elev^2, slope2 = slope^2, swhc2 = swhc^2,
                              pH2 = pH^2, expoNS2 = expoNS^2, expoEW2 = expoEW^2)
#

# function to format data for each species
# keep only the relevant SI and remove all NA
spForm <- function(df, spCode){
  temp <- df %>% select(paste0('potentiel_', spCode), elev, slope, swhc, pH,
                                         expoNS, expoEW, GRECO,
                                         elev2, slope2, swhc2, pH2, expoNS2, expoEW2) %>%
                        mutate(miss = rowSums(is.na(df))) %>%
                        filter(miss == 0) %>% select(-miss)
  return(temp)
}

# function to select model with most significant term and with AIC not higher than
# best model AIC+2
bestMod <- function(allMod){
  # first retrieve best AIC
  bestic <- summary(allMod)$bestic
  # then retrieve number of significant term in each model
  df <- data.frame()
  rank <- 1
  modAIC <- 0
  while(modAIC <= (bestic + 2)){
    res <- c(rank, nrow(summary(allMod@objects[[rank]])$coef), sum(summary(allMod@objects[[rank]])$coef[, 'Pr(>|t|)'] <= 0.05))
    df <- rbind(df, res)
    rank <- rank + 1
    modAIC <- summary(allMod@objects[[rank]])$aic
  }
  colnames(df) <- c('rank', 'nCoef', 'nCoefSign')
  df$propCoefSign <- df$nCoefSign / df$nCoef
  df <- df[order(-df$propCoefSign, df$rank), ]
  return(df)
}


###############################################################
# model sp 03
###############################################################

# data
salemSI03 <- spForm(df = salemSI, spCode = '03')

# run all possible models
# start <- Sys.time()
allMod03 <- glmulti(potentiel_03 ~ ., data = salemSI03, level = 1, method = 'h',
               crit = aic, plotty = TRUE, report = TRUE,
               family = Gamma(link = 'log'))
# end <- Sys.time()
# end - start

# plot(allMod03, type = 'p')
# plot(allMod03, type = 'w')
# plot(allMod03, type = 's')

# overview
print(allMod03)

# check best model
summary(allMod03@objects[[1]])

# select model with most significant term and with AIC not higher than
# best model AIC+2
df03 <- bestMod(allMod03)
summary(allMod03@objects[[df03[1, 'rank']]])
modQpetraea <- allMod03@objects[[df03[1, 'rank']]]
saveRDS(modQpetraea , './data/salemSI/modQpetraea.rds')

###############################################################
# model sp 09
###############################################################

# data
salemSI09 <- spForm(df = salemSI, spCode = '09')

# run all possible models
allMod09 <- glmulti(potentiel_09 ~ ., data = salemSI09, level = 1, method = 'h',
               crit = aic, plotty = TRUE, report = TRUE,
               family = Gamma(link = 'log'))
#
# check best model
summary(allMod09@objects[[1]])

# select model with most significant term and with AIC not higher than
# best model AIC+2
df09 <- bestMod(allMod09)
summary(allMod09@objects[[df09[1, 'rank']]])
modFsylvatica <- allMod09@objects[[df09[1, 'rank']]]
saveRDS(modFsylvatica , './data/salemSI/modFsylvatica.rds')

###############################################################
# model sp 61
###############################################################

# data
salemSI61 <- spForm(df = salemSI, spCode = '61')

# run all possible models
allMod61 <- glmulti(potentiel_61 ~ ., data = salemSI61, level = 1, method = 'h',
               crit = aic, plotty = TRUE, report = TRUE,
               family = Gamma(link = 'log'))
#
# check best model
summary(allMod61@objects[[1]])

# select model with most significant term and with AIC not higher than
# best model AIC+2
df61 <- bestMod(allMod61)
summary(allMod61@objects[[df61[1, 'rank']]])
modAalba <- allMod61@objects[[df61[1, 'rank']]]
saveRDS(modAalba , './data/salemSI/modAalba.rds')

###############################################################
# model sp 62
###############################################################

# data
salemSI62 <- spForm(df = salemSI, spCode = '62')

# run all possible models
allMod62 <- glmulti(potentiel_62 ~ ., data = salemSI62, level = 1, method = 'h',
               crit = aic, plotty = TRUE, report = TRUE,
               family = Gamma(link = 'log'))
#
# check best model
summary(allMod62@objects[[1]])

# select model with most significant term and with AIC not higher than
# best model AIC+2
df62 <- bestMod(allMod62)
summary(allMod62@objects[[df62[1, 'rank']]])
modPabies <- allMod62@objects[[df62[1, 'rank']]]
saveRDS(modPabies , './data/salemSI/modPabies.rds')

###############################################################
# manually add linear effects when quadratic term are present
# and remove terms becoming non-signifant because of this operation
###############################################################

# 03
summary(modQpetraea, correlation = TRUE)
# add missing linear effects and remove ns effects
modQpetraea <- glm(potentiel_03 ~ slope + swhc + expoNS + slope2 + elev + pH, family = Gamma(link = 'log'), data = salemSI03)
# save
saveRDS(modQpetraea , './data/salemSI/modQpetraea.rds')

# 09
summary(modFsylvatica, correlation = TRUE)
# add missing linear effects and remove ns effects
modFsylvatica <- (glm(potentiel_09 ~ slope + swhc + pH + slope2 + swhc2 +
pH2 + elev, family = Gamma(link = 'log'), data = salemSI09))
# save
saveRDS(modFsylvatica , './data/salemSI/modFsylvatica.rds')

# 61
summary(modAalba, correlation = TRUE)
modAalba <- glm(potentiel_61 ~ GRECO + slope + swhc + pH + expoNS + expoEW + slope2 +
swhc2 + pH2 + elev, family = Gamma(link = 'log'), data = salemSI61)
# save
saveRDS(modAalba , './data/salemSI/modAalba.rds')

# 62
summary(modPabies, correlation = TRUE)
modPabies <- glm(potentiel_62 ~ GRECO + elev + slope + swhc + expoNS + slope2 +
swhc2 + pH, family = Gamma(link = 'log'), data = salemSI62)
# save
saveRDS(modPabies , './data/salemSI/modPabies.rds')

###############################################################
# evaluation
###############################################################

# prediction
salemSI03$pred <- predict(modQpetraea, newdata = salemSI03, type = 'response')
salemSI09$pred <- predict(modFsylvatica, newdata = salemSI09, type = 'response')
salemSI61$pred <- predict(modAalba, newdata = salemSI61, type = 'response')
salemSI62$pred <- predict(modPabies, newdata = salemSI62, type = 'response')

# stack sp df
salemSI03 <- salemSI03 %>% rename(potentiel = potentiel_03) %>% mutate(sp = as.factor('03'))
salemSI09 <- salemSI09 %>% rename(potentiel = potentiel_09) %>% mutate(sp = as.factor('09'))
salemSI61 <- salemSI61 %>% rename(potentiel = potentiel_61) %>% mutate(sp = as.factor('61'))
salemSI62 <- salemSI62 %>% rename(potentiel = potentiel_62) %>% mutate(sp = as.factor('62'))
df <- rbind(salemSI03, salemSI09, salemSI61, salemSI62)

# define min and max limits for plot extent
df <- df %>% group_by(sp) %>%
             mutate(maxLim = max(potentiel, pred), minLim = min(potentiel, pred)) %>%
             ungroup()
#

# obs vs pred models
r03 <- summary(lm(potentiel ~ pred, data = df[df$sp == '03',]))$r.squared
r09 <- summary(lm(potentiel ~ pred, data = df[df$sp == '09',]))$r.squared
r61 <- summary(lm(potentiel ~ pred, data = df[df$sp == '61',]))$r.squared
r62 <- summary(lm(potentiel ~ pred, data = df[df$sp == '62',]))$r.squared
rdf <- data.frame(sp = c('03', '09', '61', '62'), r = c(r03, r09, r61, r62))

# plot
pl1 <- ggplot(data = df,aes(x = pred, y = potentiel)) +
    geom_point( ) +
    geom_abline(slope = 1, intercept = 0, lwd = 1) +
    geom_smooth(method = 'lm', formula = y ~ x) +
    geom_point(aes(x = maxLim, y = minLim), alpha = 0) +
    geom_point(aes(x = minLim, y = maxLim), alpha = 0) +
    geom_text(data = rdf, aes(x = -Inf, y = Inf, label = paste('r2=',round(r, 3))),
              hjust = -1, vjust = 3, show.legend = FALSE, col = 'red') +
    facet_wrap(. ~ sp, scale = 'free') +
    ylab('observations') +
    xlab('predictions') +
    theme_light() +
    theme(aspect.ratio = 1, panel.grid.minor = element_blank(),
                  strip.background = element_blank(),
                  strip.text = element_text(colour = 'black'),
                  legend.position = 'bottom',
                  legend.title = element_blank(),
                  panel.spacing = unit(20, 'pt'))
#
# save plot
ggsave(file = './initialLandscape/evaluation/salemSI.pdf', plot = pl1, width = 10, height = 10)


# plot residuals against variables and prediction ranges
# define prediction range of variables
predRange <- data.frame(var = c('elev', 'slope', 'swhc', 'pH', 'expoNS', 'expoEW'),
                        min = c(237, 0, 1, 4.3, -1, -1),
                        max = c(2794, 80, 14, 7, 1, 1))
#

# function to retrieve names of continuous simple effects
simpleF <- function(mod){
  simple <- names(coef(mod)[!(substr(names(coef(mod)), 1, 2) %in% c('(I', 'GR')) &
    substr(names(coef(mod)), nchar(names(coef(mod))),
    nchar(names(coef(mod)))) != '2'])
  return(simple)
}
# function to retrieve names of continuous quadratic effects
quadraF <- function(mod){
  quadra <- names(coef(mod)[substr(names(coef(mod)), nchar(names(coef(mod))),
    nchar(names(coef(mod)))) == '2'])
  return(quadra)
}
# function to retrieve names of effect present as simple AND quadratic
bothF <- function(simple, quadra){
  quadra <- substr(quadra, 1, nchar(quadra)-1)
  all <- c(simple, quadra)
  both <- all[duplicated(all)]
  return(both)
}
# function to retrieve names of effect present as simple OR quadratic
onlyF <- function(both, simple, quadra){
  quadra <- substr(quadra, 1, nchar(quadra)-1)
  all <- c(simple, quadra)
  only <- all[!(all %in% both)]
  return(only)
}

# function to plot residuals against variables and to evaluate the effetc of
# predicting out of range of calibration
plotMod <- function(name, mod, predRange, df, type = c('s', 'q', 'sq')){
  # coef(mod)[name]
  if(type == 's'){
    curve(coef(mod)[name] * x, min(df[, name], predRange[predRange$var == name, 'min']), max(df[, name], predRange[predRange$var == name, 'max']), main = name)
  } else if(type == 'q'){
    curve(coef(mod)[name] * x*x, min(df[, name], predRange[predRange$var == substr(name, 1, nchar(name)-1), 'min']), max(df[, name], predRange[predRange$var == substr(name, 1, nchar(name)-1), 'max']), main = name)
  } else if(type == 'sq'){
    curve(coef(mod)[name] * x + coef(mod)[paste0(name, 2)] *x*x, min(df[, name], predRange[predRange$var == name, 'min']), max(df[, name], predRange[predRange$var == name, 'max']), main = name)
  }
  abline(v = range(df[, name]), lty = 2)
  plot(resid(mod) ~ df[, name], main = name)
  panel.smooth(df[, name], resid(mod))
  abline(h = 0, lty = 2)
  return()
}



# function to run all at once
plotSp <- function(mod, plotMod, predRange, df){
  simple <- simpleF(mod)
  quadra <- quadraF(mod)
  both <- bothF(simple, quadra)
  only <- onlyF(both, simple, quadra)
  par(mfrow = c(length(c(only, both)), 2))
  # only simple effects
  lapply(c(only, simple)[duplicated(c(only, simple))], plotMod, mod = mod, predRange = predRange, df = df, type = 's')
  # simple and quadratic effect
  lapply(both, plotMod, mod = mod, predRange = predRange, df = df, type = 'sq')
  # only quadratic effect
  lapply(quadra[!quadra %in% paste0(both,2)], plotMod, mod = mod, predRange = predRange, df = df, type = 'q')
}

# 03
pdf('./initialLandscape/evaluation/modQpet.pdf', width = 7, height = 15)
plotSp(modQpetraea, plotMod, predRange, salemSI03)
dev.off()

# 09
pdf('./initialLandscape/evaluation/modFsyl.pdf', width = 7, height = 15)
plotSp(modFsylvatica, plotMod, predRange, salemSI09)
dev.off()

# 61
pdf('./initialLandscape/evaluation/modAalb.pdf', width = 7, height = 15)
plotSp(modAalba, plotMod, predRange, salemSI61)
dev.off()

# 62
pdf('./initialLandscape/evaluation/modPabi.pdf', width = 7, height = 15)
plotSp(modPabies, plotMod, predRange, salemSI62)
dev.off()
