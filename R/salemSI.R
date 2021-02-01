

# predict salem SI (site index) from PROTEST models
salemSI <- function(envdf){

  require(dplyr)
  # load SI models
  Qpet <- readRDS('./data/salemSImodels/modQpetraea.rds')
  Fsyl <- readRDS('./data/salemSImodels/modFsylvatica.rds')
  Aalb <- readRDS('./data/salemSImodels/modAalba.rds')
  Pabi <- readRDS('./data/salemSImodels/modPabies.rds')

  # change variable names
  envdf <- envdf %>% rename(rum = swhc, alti = elev, greco = GRECO)

  # create ExpoNS variable
  envdf$expoNS <- cos(envdf$aspect*pi/180)

  # predict
  envdf$SIQpet <- predict(Qpet, newdata = envdf, type = 'response')
  envdf$SIFsyl <- predict(Fsyl, newdata = envdf, type = 'response')
  envdf$SIAalb <- predict(Aalb, newdata = envdf, type = 'response')
  envdf$SIPabi <- predict(Pabi, newdata = envdf, type = 'response')

  # change back variable names
  envdf <- envdf %>% rename(swhc = rum, elev = alti, GRECO = greco)

  # remove ExpoNS variable
  envdf <- envdf %>% select(-expoNS)

  return(envdf)

}
