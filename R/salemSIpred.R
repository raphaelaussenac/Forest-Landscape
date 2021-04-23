# prediction algo
#
# In the calibration data slope ranges from 0 to 47 while it ranges from 0 to 80
# in the prediction data. Besides in the SI models, slope has a quadratic effect.
# Predictions out of the calibration range might therefore not be reliable. We
# therefore set the max slope in the prediction data to the max slope in the
# calibration data.

# predict salem SI (site index) from PROTEST models
salemSI <- function(df){
  require(dplyr)

  dfBackup <- df

  # set maximum slope to max slope of calibration data
  df[!is.na(df$slope) & df$slope > 47.72631, 'slope'] <- 47.72631

  # load SI models
  Qpet <- readRDS('./data/salemSI/modQpetraea.rds')
  Fsyl <- readRDS('./data/salemSI/modFsylvatica.rds')
  Aalb <- readRDS('./data/salemSI/modAalba.rds')
  Pabi <- readRDS('./data/salemSI/modPabies.rds')

  # manage variable class
  df$forest <- as.factor(df$forest)
  df$GRECO <- as.factor(df$GRECO)

  # create missing variables
  df <- df %>% mutate(expoNS = cos(df$aspect*pi/180), expoEW = sin(df$aspect*pi/180),
                      elev2 = elev^2, slope2 = slope^2, swhc2 = swhc^2,
                      pH2 = pH^2, expoNS2 = expoNS^2, expoEW2 = expoEW^2)
  #

  # predict SI only in forest cells
  df <- df %>% filter(forest == 1)
  df$SIQpet <- predict(Qpet, newdata = df, type = 'response')
  df$SIFsyl <- predict(Fsyl, newdata = df, type = 'response')
  df$SIAalb <- predict(Aalb, newdata = df, type = 'response')
  df$SIPabi <- predict(Pabi, newdata = df, type = 'response')

  # add SI to df
  df <- merge(dfBackup, df[, c('cellID25', 'SIQpet', 'SIFsyl', 'SIAalb', 'SIPabi')],
              by = 'cellID25', all.x = TRUE)
  #

  return(df)

}
