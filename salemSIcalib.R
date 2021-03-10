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
setwd("C:/Users/raphael.aussenac/Documents/GitHub/LandscapeInit")

# load species-specific SI and env variables at NFI plots
salemSI <- read.csv('./data/salemSI/bdBauges_for_SI_calibration_2021_03_09.txt', sep = '\t')

# load park limits
park <- readOGR(dsn = "./data/GEO", layer = "park", encoding = "UTF-8", use_iconv = TRUE)

# load geological and pH data
Cd_crbn <- raster("./data/init/Cd_crbn.asc")
Cd_hydr <- raster("./data/init/Cd_hydr.asc")
pH <- raster("./initialLandscape/pH.asc")

###############################################################
# assign geol and pH values to all NFI plots
###############################################################

# convert NFI into spatial points
NFIplots <- SpatialPointsDataFrame(salemSI[,c("xl93", "yl93")], data = data.frame(salemSI[,'idp']), proj4string = CRS(proj4string(park)))

# extract geol and pH values at NFI plots
salemSI$Cd_crbn <- extract(Cd_crbn, NFIplots)
salemSI$Cd_hydr <- extract(Cd_hydr, NFIplots)
salemSI$pH <- extract(pH, NFIplots)

###############################################################
# manage data format, variable class, unit, model selection...
###############################################################

# remove useless variables
salemSI <- salemSI %>% select(-idp, -xl93, -yl93, -aspect)

# define variable class
salemSI$GRECO <- as.factor(salemSI$GRECO)
salemSI$Cd_crbn <- as.factor(salemSI$Cd_crbn)
salemSI$Cd_hydr <- as.factor(salemSI$Cd_hydr)

# convert slope from degrees to % because SI where built with % and
# the relationship between degrees and percent is not linear
salemSI$slope <- tan(salemSI$slope * (2*pi) /360) * 100

# add quadratic effects
salemSI <- salemSI %>% mutate(elev2 = elev^2, slope2 = slope^2, swhc2 = swhc^2,
                              pH2 = pH^2, expoNS2 = expoNS^2, expoEW2 = expoEW^2)
#

# function to format data for each species
# keep only the relevant SI and remove all NA
spForm <- function(df, spCode){
  temp <- df %>% select(paste0('potentiel_', spCode), elev, slope, swhc, pH,
                                         expoNS, expoEW, GRECO, Cd_crbn, Cd_hydr,
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
    res <- c(rank, nrow(summary(allMod@objects[[rank]])$coef), sum(summary(allMod@objects[[rank]])$coef[, 'Pr(>|t|)'] <= 0.1))
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
               family = Gamma(link = "log"))
# end <- Sys.time()
# end - start

# plot(allMod03, type = "p")
# plot(allMod03, type = "w")
# plot(allMod03, type = "s")

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
               family = Gamma(link = "log"))
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
               family = Gamma(link = "log"))
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
               family = Gamma(link = "log"))
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
    facet_wrap(. ~ sp, scale = "free") +
    ylab('observations') +
    xlab('predictions') +
    theme_light() +
    theme(aspect.ratio = 1, panel.grid.minor = element_blank(),
                  strip.background = element_blank(),
                  strip.text = element_text(colour = 'black'),
                  legend.position = "bottom",
                  legend.title = element_blank(),
                  panel.spacing = unit(20, 'pt'))
#
# save plot
ggsave(file = './initialLandscape/evaluation/salemSI.pdf', plot = pl1, width = 10, height = 10)
