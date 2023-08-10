# clean up environment
rm(list = ls())

# set work directory
setwd('~/Documents/code/Forest-Landscape')

# load packages
library(ggplot2)
library(dplyr)
library(tidyr)

# select landscape (bauges, milicz, sneznik)
landscape <- 'sneznik'

# open new tree data and published tree data
new <- read.csv(paste0('./', landscape, '/trees.csv'))
old <- read.csv(paste0('./I-MAESTRO_data/', landscape, '/trees.csv'))
# TODO 1: compare with published data

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
# plot
ggplot(data = spDiag) +
geom_histogram(aes(x = diff_pct)) +
theme_bw()
# list of species with variables deviation >= X%
X = 1
strongDeviat <- spDiag %>% filter(abs(diff_pct) >= X )
strongDeviat
spList <- unique(strongDeviat$sp)
spList

########################################################################
# Plot species frequency and show deviating sp
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

########################################################################
# Plot species dbh and h deviations
# define new df
dfDeviat <- df %>% filter(sp %in% spList)

# dbh distribution
ggplot(data = dfDeviat) +
geom_density(aes(x = dbh, fill = dataset), alpha = 0.5) +
facet_grid(sp ~ ., scales = 'free') +
theme_bw()

# h distribution
ggplot(data = dfDeviat) +
geom_density(aes(x = h, fill = dataset), alpha = 0.5) +
facet_grid(sp ~ ., scales = 'free') +
theme_bw()





# TODO: diff old vs new of all output rasters
# ou comparaison cellID?