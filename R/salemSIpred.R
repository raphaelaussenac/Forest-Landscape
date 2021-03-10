# predict salem SI (site index) from PROTEST models
salemSI <- function(df){

  dfBackup <- df

  require(dplyr)
  # load SI models
  Qpet <- readRDS('./data/salemSI/modQpetraea.rds')
  Fsyl <- readRDS('./data/salemSI/modFsylvatica.rds')
  Aalb <- readRDS('./data/salemSI/modAalba.rds')
  Pabi <- readRDS('./data/salemSI/modPabies.rds')

  # manage variable class
  df$forest <- as.factor(df$forest)
  df$GRECO <- as.factor(df$GRECO)
  df$Cd_crbn <- as.factor(df$Cd_crbn)
  df$Cd_hydr <- as.factor(df$Cd_hydr)

  # make sure prediction are made on same unit than calibration
  # convert slope degrees in %
  df$slope <- tan(df$slope * (2*pi) /360) * 100

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
  df <- merge(dfBackup, df[, c('cellID', 'SIQpet', 'SIFsyl', 'SIAalb', 'SIPabi')],
              by = 'cellID', all.x = TRUE)
  #

  return(df)

}
