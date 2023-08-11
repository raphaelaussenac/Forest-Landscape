# clean up environment
rm(list = ls())

# set work directory
setwd('~/Documents/code/Forest-Landscape')

# load packages
library(ggplot2)
library(dtplyr)
library(dplyr)
library(tidyr)
library(terra)

# select landscape (bauges, milicz, sneznik)
landscape <- 'bauges'

# open new tree data and published tree data
new <- read.csv(paste0('./', landscape, '/trees.csv'))
old <- read.csv(paste0('./I-MAESTRO_data/', landscape, '/trees.csv'))

# stack data
new$dataset <- 'new'
old$dataset <- 'old'
df <- rbind(new, old)

########################################################################
# Diagnostics function
########################################################################

diagnostics <- function(level = c('dataset', 'sp'), df){
        if(level == 'dataset'){
                diag <- df %>% group_by(dataset)
        }
        if(level == 'sp'){
                diag <- df %>% group_by(dataset, sp)
        }
        diag <- diag %>%
               summarise(nsp = length(unique(sp)),
                         ntot = sum(n),
                         Qd0 = quantile(dbh, probs = 0),
                         Qd25 = quantile(dbh, probs = 0.25),
                         Qd50 = quantile(dbh, probs = 0.5),
                         Qd75 = quantile(dbh, probs = 0.75),
                         Qd100 = quantile(dbh, probs = 1),
                         Qh0 = quantile(h, probs = 0),
                         Qh25 = quantile(h, probs = 0.25),
                         Qh50 = quantile(h, probs = 0.5),
                         Qh75 = quantile(h, probs = 0.75),
                         Qh100 = quantile(h, probs = 1))
        if(level == 'dataset'){
                i <- 2
        }
        if(level == 'sp'){
                i <- 3
        }
        diag <- diag %>% pivot_longer(cols = colnames(diag)[i:length(colnames(diag))])
        diag <- diag %>% pivot_wider(names_from = dataset)
        diag <- diag %>% mutate(diff_pct = round((new * 100 / old) - 100, 2))
        return(diag)
}


########################################################################
# General Diagnostics
genDiag <- diagnostics('dataset', df)
genDiag

########################################################################
# Species-level Diagnostics
spDiag <- diagnostics('sp', df) %>% arrange(-abs(diff_pct))
spDiag

# number of individuals per species
spDiag_ntot <- spDiag %>% filter(name == 'ntot')

# plot
ggplot(data = spDiag_ntot) +
geom_histogram(aes(x = diff_pct)) +
theme_bw()
# list of species with ntot deviation >= X%
X = 1
strongDeviat <- spDiag_ntot %>% filter(abs(diff_pct) >= X )
strongDeviat
spList <- unique(strongDeviat$sp)
spList

# Plot species frequency and show deviating sp
if(length(spList) > 0){

        freq <- df %>% group_by(dataset, sp) %>%
               summarise(n = sum(n)) %>%
               mutate(deviatSp = as.factor(case_when(sp %in% spList ~ 'deviation', TRUE ~ 'ok')))

        #
        ggplot(data = freq) +
        geom_bar(aes(x = reorder(sp, -n), y = n, fill = dataset, col = deviatSp), alpha = 0.1, width = 0.5, stat = 'identity', position = position_dodge()) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
                plot.title = element_text(hjust = 0.5)) +
        ggtitle('sp deviating in red')
}

########################################################################
# Plot dbh distribution of deviating species 
# list of species with diameter quantile deviation >= X%
X = 1
strongDeviat <- spDiag %>% filter(substr(spDiag$name,1,2) == 'Qd', abs(diff_pct) >= X )
strongDeviat
spList <- unique(strongDeviat$sp)
spList

# define new df
dfDeviat <- df %>% filter(sp %in% spList)

# dbh distribution
ggplot(data = dfDeviat) +
geom_density(aes(x = dbh, fill = dataset), alpha = 0.5) +
facet_grid(sp ~ ., scales = 'free') +
theme_bw()

########################################################################
# Plot h distribution of deviating species 
# list of species with diameter quantile deviation >= X%
X = 1
strongDeviat <- spDiag %>% filter(substr(spDiag$name,1,2) == 'Qh', abs(diff_pct) >= X )
strongDeviat
spList <- unique(strongDeviat$sp)
spList

# define new df
dfDeviat <- df %>% filter(sp %in% spList)
# h distribution
ggplot(data = dfDeviat) +
geom_density(aes(x = h, fill = dataset), alpha = 0.5) +
facet_grid(sp ~ ., scales = 'free') +
theme_bw()

########################################################################
# spatial diag
########################################################################

# open cellID raster
cellID25 <- rast(paste0('./', landscape, '/cellID25.asc'))
celldf <- as.data.frame(cellID25)

# summarise values by cell
df <- lazy_dt(df)
syn <- df %>% group_by(dataset, cellID25) %>%
               summarise(BA = sum((pi * (dbh/200)^2) * n),
                         meanD = mean(dbh),
                         meanH = mean(h),
                         n = sum(n),
                         nsp = length(unique(sp))) %>%
                         as.data.frame()
# merge data
synNew <- left_join(celldf, syn %>% filter(dataset == 'new'), by = 'cellID25')
synOld <- left_join(celldf, syn %>% filter(dataset == 'old'), by = 'cellID25')

# BA
cellID25$BAnew <- synNew$BA
cellID25$BAold <- synOld$BA

# meanD
cellID25$meanDnew <- synNew$meanD
cellID25$meanDold <- synOld$meanD

# meanH
cellID25$meanHnew <- synNew$meanH
cellID25$meanHold <- synOld$meanH

# n
cellID25$nnew <- synNew$n
cellID25$nold <- synOld$n

# nsp
cellID25$nspnew <- synNew$nsp
cellID25$nspold <- synOld$nsp

# diff
cellID25$BAdiff <- cellID25$BAnew - cellID25$BAold
cellID25$meanDdiff <- cellID25$meanDnew - cellID25$meanDold
cellID25$meanHdiff <- cellID25$meanHnew - cellID25$meanHold
cellID25$ndiff <- cellID25$nnew - cellID25$nold
cellID25$nspdiff <- cellID25$nspnew - cellID25$nspold

# plot
diff <- cellID25[[12:16]]
plot(diff)