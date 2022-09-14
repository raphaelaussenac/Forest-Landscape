###############################################################
# initialisation
###############################################################

# clean up environment
rm(list = ls())

# set work directory
setwd('~/Documents/code/Forest-Landscape')

# load packages
library(dplyr)
library(ggplot2)
library(gridExtra)
library(stringr)
library(ggridges)
library(hrbrthemes)
library(tidyr)

# load tree data
bauges <- read.csv('./bauges/trees.csv') %>% mutate(site = 'bauges')
milicz <- read.csv('./milicz/trees.csv') %>% mutate(site = 'milicz')
sneznik <- read.csv('./sneznik/trees.csv') %>% mutate(site = 'sneznik')
df <- bind_rows(bauges, milicz, sneznik)

# create synthesis folder
if (!(dir.exists('./synthesis'))) {dir.create('./synthesis', recursive = TRUE)}

###############################################################
# synthesis table
###############################################################

syndf <- df %>% group_by(site, sp) %>%
            summarise(n = sum(n),
                      dbhQ5 = quantile(dbh, probs = 0.05),
                      dbhQ25 = quantile(dbh, probs = 0.25),
                      dbhQ50 = quantile(dbh, probs = 0.5),
                      dbhQ75 = quantile(dbh, probs = 0.75),
                      dbhQ95 = quantile(dbh, probs = 0.95),
                      hQ5 = quantile(h, probs = 0.05),
                      hQ25 = quantile(h, probs = 0.25),
                      hQ50 = quantile(h, probs = 0.5),
                      hQ75 = quantile(h, probs = 0.75),
                      hQ95 = quantile(h, probs = 0.95),
                      BAsp = sum((pi * (dbh/200)^2) * n)) %>%
            group_by(site) %>%
            mutate(BAtot = sum(BAsp),
                   BAprop = round(BAsp*100/BAtot, 4)) %>%
            select(-BAsp, -BAtot) %>%
            relocate(BAprop, .after = n) %>%
            group_by(site, sp) %>%
            arrange(-BAprop) %>%
            ungroup()

write.csv(syndf %>% filter(site == 'bauges'), './synthesis/bauges.csv', row.names = FALSE)
write.csv(syndf %>% filter(site == 'milicz'), './synthesis/milicz.csv', row.names = FALSE)
write.csv(syndf %>% filter(site == 'sneznik'), './synthesis/sneznik.csv', row.names = FALSE)


###############################################################
# main sp donut G
###############################################################

mainG <- df %>% group_by(site, sp) %>%
                 summarise(BAsp = sum((pi * (dbh/200)^2) * n)) %>%
                 group_by(site) %>%
                 mutate(BAtot = sum(BAsp),
                        prop = BAsp*100/BAtot,
                        sp = ifelse(prop>=5, sp, 'other')) %>%
                 group_by(site, sp) %>%
                 summarise(prop = sum(prop)) %>%
                 arrange(site, -prop) %>%
                 group_by(site) %>%
                 mutate(ymax = cumsum(prop),
                        ymin = c(0, head(ymax, n=-1))) %>%
                 group_by(site, sp)
mainG <- mainG %>% mutate(sp1 = factor(sp, levels = sort(unique(mainG$sp), decreasing = TRUE)),
                        site  = str_to_title(site))

mainG2 <- mainG %>% mutate(label = (ymax - ymin) / 2,
                              sp = str_replace(sp, ' ', ' \n '))

# main species donut
pl1 <- ggplot(data = mainG2, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = sp1)) +
geom_rect() +
coord_polar(theta = 'y') +
xlim(c(2, 5)) +
theme_void() +
theme(legend.position = 'none', legend.title = element_blank(),
      strip.text.x = element_blank(), panel.spacing.y = unit(-3, "lines")) + #panel.spacing.y = unit(-3, "lines")
facet_wrap( ~ site, ncol = 1) +
scale_fill_brewer(palette = "Spectral")

# species labels
pl1 <- pl1 + geom_text(
    aes(label = sp,
    x = 5,
    y = (ymin + ymax) / 2),
    # check_overlap = TRUE,
    size = 5,
    inherit.aes = TRUE,
    show.legend = FALSE
  )

# proportion labels
pl1 <- pl1 + geom_text(
    aes(label = paste0(round(prop,1), '%'),
    x = 3.5,
    y = (ymin + ymax) / 2),
    # check_overlap = TRUE,
    size = 4,
    inherit.aes = TRUE,
    show.legend = FALSE
  )

# central title
pl1 <-  pl1 + geom_text(
    aes(label = str_to_title(site),
    x = rep(2, nrow(mainG2)),
    y = rep(2, nrow(mainG2))),
    # check_overlap = TRUE,
    size = 7,
    inherit.aes = TRUE,
    show.legend = FALSE
  )

# arrows
pl1 <- pl1 + geom_segment(aes(x = rep(4, nrow(mainG2)), y = (ymin + ymax) / 2, xend =rep(4.3, nrow(mainG2)), yend = (ymin + ymax) / 2), inherit.aes = TRUE)
pl1 <- pl1 + ggtitle('proportion of species total basal area') + theme(plot.title = element_text(hjust = 0.5, size = 17, face = "bold"))
ggsave(file = './synthesis/spG.pdf', plot = pl1, width = 15, height = 12)

###############################################################
# main sp dbh and h distribution
###############################################################

mainsp <- df %>% filter(sp %in% mainG$sp)
# set main sp to others deending on sites
mainsp[mainsp$site == 'bauges' & mainsp$sp %in% c('Pinus sylvestris', 'Quercus robur'), 'sp'] <- 'other'
mainsp[mainsp$site == 'milicz' & mainsp$sp %in% c('Abies alba', 'Fraxinus excelsior', 'Picea abies'), 'sp'] <- 'other'
mainsp[mainsp$site == 'sneznik' & mainsp$sp %in% c('Pinus sylvestris'), 'sp'] <- 'other'
raresp <- df %>% filter(!(sp %in% mainG$sp)) %>% mutate(sp = 'other')
dist <- bind_rows(mainsp, raresp)
dist <- dist %>% mutate(sp = factor(sp, levels = sort(unique(dist$sp), decreasing = TRUE)),
                    site  = str_to_title(site))
dist2 <- dist[rep(seq_along(dist$n), dist$n), ] # weigtht dbh distributions by n
dist2 <- pivot_longer(data = dist2,
                     cols = c(dbh,h))

dist3 <- dist2 %>% group_by(site, sp, name) %>% summarise(dbh = mean(value))

# plot dbh
pl2 <- ggplot(dist2 %>% filter(name == 'dbh'), aes(x = value, y = sp, fill = sp, col = sp)) +
stat_density_ridges() +
scale_fill_brewer(palette = "Spectral") +
scale_colour_brewer(palette = "Spectral") +
theme_bw() +
theme(legend.position = "none") +
geom_text(data = dist3,
    aes(label = sp,
    x = 70,
    y = sp),
    col = 'black',
    # check_overlap = TRUE,
    size = 5, vjust = -1, hjust = 1,
    inherit.aes = TRUE,
    show.legend = FALSE
  ) +
facet_wrap(~ site, scales = 'free_y', ncol = 1) +
scale_y_discrete(position = "right") +
xlim(0, 80) +
theme(strip.background = element_blank(), strip.text = element_blank(), axis.title = element_blank(),
        axis.ticks.y = element_blank(), panel.border = element_blank(), axis.text.y = element_blank()) #, panel.spacing.y = unit(3, "lines")
pl2 <- pl2 + ggtitle('distribution of species diameter (cm)') + theme(plot.title = element_text(hjust = 0.5, size = 17, face = "bold"))
ggsave(file = './synthesis/dbhmainsp.pdf', plot = pl2, width = 15, height = 12)


# plot h
pl3 <- ggplot(dist2 %>% filter(name == 'h'), aes(x = value, y = sp, fill = sp, col = sp)) +
stat_density_ridges() +
scale_fill_brewer(palette = "Spectral") +
scale_colour_brewer(palette = "Spectral") +
theme_bw() +
theme(legend.position = "none") +
geom_text(data = dist3,
    aes(label = sp,
    x = 40,
    y = sp),
    col = 'black',
    # check_overlap = TRUE,
    size = 5, vjust = -1, hjust = 1,
    inherit.aes = TRUE,
    show.legend = FALSE
  ) +
facet_wrap(~ site, scales = 'free_y', ncol = 1) +
scale_y_discrete(position = "right") +
xlim(0, 45) +
theme(strip.background = element_blank(), strip.text = element_blank(), axis.title = element_blank(),
        axis.ticks.y = element_blank(), panel.border = element_blank(), axis.text.y = element_blank()) #, panel.spacing.y = unit(3, "lines")
pl3 <- pl3 + ggtitle('distribution of species height (m)') + theme(plot.title = element_text(hjust = 0.5, size = 17, face = "bold"))
ggsave(file = './synthesis/hmainsp.pdf', plot = pl3, width = 10, height = 10)

# assemble plot
pl4 <- grid.arrange(pl1, pl2, pl3, ncol = 3, nrow = 1, widths = c(0.25, 0.3, 0.3)) #, heights = c(1, 2))
ggsave(file = './synthesis/synth.pdf', plot = pl4, width = 15, height = 12)
