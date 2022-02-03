###############################################################
# initialisation
###############################################################

# clean up environment
rm(list = ls())

# set work directory
setwd('./Documents/github/Forest-Landscape/temp')
# TODO: no wd --> directly from evalPath

# load library
library(ggplot2)
library(dplyr)

# load msd data
bauges <- read.csv('msd_bauges.csv')
milicz <- read.csv('msd_milicz.csv')
sneznik <- read.csv('msd_sneznik.csv')

# merge data
df <- bind_rows(bauges, milicz, sneznik) %>% filter(var != 'msd')
df$landscape <- factor(df$landscape, levels = c('bauges', 'sneznik', 'milicz'))

# plot
pl1 <- ggplot(df) +
geom_bar(aes(x = landscape, y = components, fill = var), stat = 'identity', position = 'stack') +
ylab('mean square deviation') +
theme_bw() +
theme(legend.position = 'bottom') +
guides(fill = guide_legend(title="MSD components"))
ggsave(file = '../heightMSD.pdf', plot = pl1, width = 5, height = 10)
