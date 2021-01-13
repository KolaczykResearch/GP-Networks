# Functions ----
source("/R/gp_functions.survival.R")
source("/R/dist.R")

# Simulation ----
set.seed(575)
m <- 100
Y1 <- rnorm(m/2, 2 , 0.8)
Y2 <- rnorm(m/2, 4, 1)
# Y1 <- rnorm(m/2, 3, 0.8)
# Y2 <- ifelse(runif(m/2) < .4, rnorm(m/2, 4, 1), rnorm(m/2, 2, .8))
Y <- abs(c(Y1, Y2))

n <- 50
G0 <- replicate(m/2, as_adj(erdos.renyi.game(n, .3, type = "gnp")), simplify = FALSE)
G1 <- replicate(m/2, as_adj(erdos.renyi.game(n, .7, type = "gnp")), simplify = FALSE)

G <- c(G0, G1)

X <- c(rep(0, m/2), rep(1, m/2))

# GP-F ----
D.F <- dist.frobenius(list(G), as.matrix(Y))
a <- c(1, 0, 5, 5)
b <- c(0, 0, 1, 1)
ns <- 1000
system.time(surv.GP_F <- gp.survival(D.F, Y, a, b, N = ns))

p1 <- plot.surv(surv.GP_F, Y, groups = X, km = F)
# p2 <- plot.surv(surv.GP_F, Y, groups = X, km = F)

# Visualize ----
gridExtra::grid.arrange(p1, p2, ncol = 2)

# add true curves in plot.surv()

# easy
# p1 <- p1 + 
#   geom_line(data = data.frame(grid
#                               , surv = 1 - pnorm(grid, 2, .8))
#             , aes(x = grid, y = surv), inherit.aes = FALSE) +
#   geom_line(data = data.frame(grid
#                               , surv = 1 - pnorm(grid, 4, 1))
#             , aes(x = grid, y = surv), inherit.aes = FALSE)

# hard
# p1 <- p1 + 
#   geom_line(data = data.frame(grid
#                               , surv = 1 - pnorm(grid, 3, .8))
#             , aes(x = grid, y = surv), inherit.aes = FALSE) +
#   geom_line(data = data.frame(grid
#                               , surv = 1 - (0.4*pnorm(grid, 4, 1) + 0.6*pnorm(grid, 2, 0.8)))
#             , aes(x = grid, y = surv), inherit.aes = FALSE)