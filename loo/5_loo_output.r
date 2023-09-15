library(gridExtra)




dbhQ <- readRDS(paste0('./loo/', landscape, '_dbhQ.RDS'))
Hdom <- readRDS(paste0('./evalHeight/heights', '_', landscape, '.rds'))
hQ5 <- readRDS(paste0('./loo/', landscape, '_0.5', '_hQ.RDS'))
hQ95 <- readRDS(paste0('./loo/', landscape, '_0.95', '_hQ.RDS'))
if (landscape == 'bauges'){
    hQ5 <- Hdom
    hQ95 <- Hdom
}
quan <- arrangeGrob(dbhQ, hQ5, hQ95, ncol = 3, nrow = 1)
stru <- readRDS(paste0('./loo/', landscape, '_loo.RDS'))
compo <- readRDS(paste0('./loo/', landscape, '_compo.RDS'))
pl <- arrangeGrob(stru, quan, compo, ncol = 1, nrow = 3)
plot(pl)

ggsave(file = paste0('./loo/', landscape, '_eval.jpg'), plot = pl, width = 20, height = 20)





# # aggregated variables
# b <- readRDS('./loo/bauges_loo.RDS')
# m <- readRDS('./loo/milicz_loo.RDS')
# s <- readRDS('./loo/sneznik_loo.RDS')

# pl <- arrangeGrob(b, m, s, ncol = 3, nrow = 1)
# ggsave(file = './loo/loo.pdf', plot = pl, width = 15, height = 15)

# # composition


# # distribution