###############################################################
# initialisation
###############################################################

# clean up environment
rm(list = ls())

# set work directory
setwd('./Documents/github/Forest-Landscape/evalHeight')

# load library
library(ggplot2)
library(dplyr)

# load msd data
bauges <- read.csv('bauges_msd.csv')
milicz <- read.csv('milicz_msd.csv')
sneznik <- read.csv('sneznik_msd.csv')

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
ggsave(file = './msd.pdf', plot = pl1, width = 5, height = 10)

# plot relative msd
dfFreq <- df %>% group_by(landscape, var) %>%
                 summarise(msd = sum(components)) %>%
                 mutate(freq = msd / sum(msd))
#
pl2 <- ggplot(dfFreq) +
geom_bar(aes(x = landscape, y = freq, fill = var), stat = 'identity', position = 'stack') +
ylab('mean square deviation') +
theme_bw() +
theme(legend.position = 'bottom') +
guides(fill = guide_legend(title="MSD components"))
