source("/R/gp_functions.gibbs_slice.R")

# Data ----

set.seed(575)
N0 <- 10; N1 <- 90
n <- 50
p0 <- .5                      # ER
pm1 <- cbind( c(.55, .45)
              , c(.45, .55) ) # SBM
G0 <- replicate(N0, as_adj(erdos.renyi.game(n, p0, type = "gnp")), simplify = FALSE)
G1 <- replicate(N1, as_adj(sample_sbm(n, pm1, c(n/2, n/2))))
G <- c(G0, G1)                        
Y <- c(rep(1, N0), rep(-1, N1))

# Unbalanced split ----
split_index <- c(rep("validate", 25), rep("train", 75)) # OCC

# GP ----
D <- kern.dis(list(G))
a <- c(10, 10)
b <- c(50, 50)

ns <- 10000
nchains <- 1
N.train <- sum(split_index == "train")
p <- dim(D)[3]
params <- list(N.train + p + 2)
for (i in 1:N.train) {
  params[i] <- paste0("f", i)
}
params[N.train + 1] <- "sigma"
for (pp in 1:p) {
  params[N.train + 1 + pp]  <- paste0("l", pp)
}
params[N.train + p + 2] <- "lp.f"
sims <- mcmc_array(ns, nchains = nchains, params)

for (ic in 1:nchains)
  sims[, ic, ] <- t(gp.class(dist = D, Y = Y, split_index = split_index
                             , a = a, b = b, ns = ns, monitor = FALSE))

# mixing
# plot(sims[, 1 ,][, N.train + 1], type = "l") # sigma
# plot(sims[, 1 ,][, N.train + 2], type = "l") # ell
# plot(sims[, 1 ,][, 1], type = "l")           # f


# OCC ----
source("/R/gp_functions.occ.R")

occ.scores <- gp.occ(sims[,1,], burn = .2, thin = 10)

# Figure 5 ----
plot.occ(occ.scores)