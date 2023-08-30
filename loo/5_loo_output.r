
# aggregated variables
b <- readRDS('./loo/bauges_loo.RDS')
m <- readRDS('./loo/milicz_loo.RDS')
s <- readRDS('./loo/sneznik_loo.RDS')

pl <- arrangeGrob(b, m, s, ncol = 3, nrow = 1)
ggsave(file = './loo/loo.pdf', plot = pl, width = 15, height = 15)

# composition


# distribution