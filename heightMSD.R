###############################################################
# initialisation
###############################################################

# clean up environment
rm(list = ls())

# set work directory
setwd('./Documents/code/Forest-Landscape/evalHeight')

# load library
library(ggplot2)
library(dplyr)
library(stringr)

# load msd data
bauges <- read.csv('bauges_msd.csv')
milicz <- read.csv('milicz_msd.csv')
sneznik <- read.csv('sneznik_msd.csv')

# merge data
df <- bind_rows(bauges, milicz, sneznik) %>% filter(var != 'msd')
df$landscape <- factor(df$landscape, levels = c('bauges', 'sneznik', 'milicz'))

# format
df$var <- toupper(df$var)
df$var <- factor(df$var, levels = c('LC', 'NU', 'SB'))
df$landscape <- str_to_title(df$landscape)
df$landscape <- factor(df$landscape, levels = c('Bauges', 'Sneznik', 'Milicz'))

# plot
pl1 <- ggplot(df) +
geom_bar(aes(x = landscape, y = components, fill = var), stat = 'identity', position = 'stack') +
scale_fill_grey(start = 0.2, end = 0.9) +
theme_bw() +
scale_y_continuous(expand = c(0, 0)) +
scale_x_discrete(expand = c(0.25, 0)) +
theme(legend.position = 'bottom', panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.ticks.x = element_blank(),
      # axis.line.x = element_line(),
      # axis.line.y = element_line()
      axis.title.y = element_blank(),
      axis.title.x = element_blank()) +
guides(fill = guide_legend(title=""))
ggsave(file = './msd.pdf', plot = pl1, width = 2, height = 4)
pl1

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
