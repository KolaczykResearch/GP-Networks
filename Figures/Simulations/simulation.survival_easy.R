# Functions ----
source("/R/gp_functions.survival.R")

# Simulation ----
set.seed(575)
n.sim <- 50
Y1 <- rnorm(n.sim, 2 , 0.8)
Y2 <- rnorm(n.sim, 4, 1)
Y <- abs(c(Y1, Y2))

G0 <- replicate(n.sim, as_adj(erdos.renyi.game(50, .3, type = "gnp")), simplify = FALSE)
G1 <- replicate(n.sim, as_adj(erdos.renyi.game(50, .7, type = "gnp")), simplify = FALSE)

G <- c(G0, G1)

D <- kern.dis(list(G), as.matrix(Y))
a <- c(1, 0, 5, 5) # flat prior on Omega and sigma
b <- c(0, 0, 1, 1)

X <- c(rep(0, n.sim), rep(1, n.sim))

N <- 1000
system.time(test <- gp.survival(D, Y, a, b, N = N))

# Check mixing
plot(test$Omega, typ = "l", ylab = expression(Omega), xlab = "iteration"
     , main = bquote(Omega ~ " ~ Exp(" ~ .(a[1]) ~ "," ~ .(b[1]) ~ ")"))
plot(test$hyperparams[1, ], typ = "l", ylab = expression(sigma), xlab = "iteration"
     , main = bquote(sigma ~ " ~ IG(" ~ .(a[2]) ~ "," ~ .(b[2]) ~ ")"))
plot(test$hyperparams[2, ], typ = "l", ylab = expression(l[net]), xlab = "iteration"
     , main = bquote(l[net] ~ " ~ IG(" ~ .(a[3]) ~ "," ~ .(b[3]) ~ ")"))
plot(test$hyperparams[3, ], typ = "l", ylab = expression(l[time]), xlab = "iteration"
     , main = bquote(l[time] ~ " ~ IG(" ~ .(a[4]) ~ "," ~ .(b[4]) ~ ")"))
plot(sapply(test$l, function(l) l[1]), typ = "l") # l(1)
plot(sapply(1:dim(test$surv)[3]
            , function(q) test$surv[1, 1, q]), typ = "l") # l(grid(1), 1)

# Figure 6 ----
plot.surv(test, Y, groups = X, km = F)

# add true curves to plot.surv
# p1 <- p1 + 
#   geom_line(data = data.frame(grid
#                               , surv = 1 - pnorm(grid, 2, .8))
#             , aes(x = grid, y = surv), inherit.aes = FALSE) +
#   geom_line(data = data.frame(grid
#                               , surv = 1 - pnorm(grid, 4, 1))
#             , aes(x = grid, y = surv), inherit.aes = FALSE)