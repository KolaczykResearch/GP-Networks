# Version : R 3.6.0
# Libraries ----
library(kernlab)
library(graphkernels)

# Functions ----
source("/R/gp_functions.gibbs_slice.R")
source("/R/dist.R")
source("/R/gp_rw.gibbs.R")
source("/R/gp_functions.occ.R")

# Data ----
set.seed(575)

m <- 100
m0 <- .1*m; m1 <- .9*m
n <- 100

pm0 <- cbind(c(.05, .15)
             , c(.15, .05))
pm1 <- cbind(c(.1, .15)
             , c(.15, .05))
G0 <- replicate(m0, as_adj(sample_sbm(n, pm0, c(n/2, n/2))))
G1 <- replicate(m1, as_adj(sample_sbm(n, pm1, c(n/2, n/2))))

G <- c(G0, G1)                        
Y <- c(rep(1, m0), rep(-1, m1))

# Unbalanced split ----
split_index <- c(rep("validate", .25*m), rep("train", .75*m)) # OCC

# GP-F ----
# GP with squared-exponential kernel and Frobenius distance
D <- dist.frobenius(list(G))
a <- c(10, 10)
b <- c(50, 50)

ns <- 2500
nchains <- 1
m.train <- sum(split_index == "train")
p <- dim(D)[3]
params <- list(m.train + p + 2)
for (i in 1:m.train) {
  params[i] <- paste0("f", i)
}
params[m.train + 1] <- "sigma"
for (pp in 1:p) {
  params[m.train + 1 + pp]  <- paste0("l", pp)
}
params[m.train + p + 2] <- "lp.f"
sims <- mcmc_array(ns, nchains = nchains, params)

for (ic in 1:nchains)
  sims[, ic, ] <- t(gp.class(dist = D, Y = Y, split_index = split_index
                             , a = a, b = b, ns = ns, monitor = FALSE))
pred.GP_F <- gp.class_pred(sims[, 1 ,], D, split_index, p
                           , burn = .2, thin = 10, avg = TRUE)
occ.GP_F <- gp.occ(sims[,1,], D, burn = .2, thin = 10)

# GP-Lambda ----
# GP with squared-exponential kernel and spectral-based distance
D <- dist.eigen(list(G), normalized = TRUE)
a <- c(10, 10)
b <- c(50, 50)

ns <- 2500
nchains <- 1
N.train <- sum(split_index == "train")
p <- dim(D)[3]
params <- list(m.train + p + 2)
for (i in 1:m.train) {
  params[i] <- paste0("f", i)
}
params[m.train + 1] <- "sigma"
for (pp in 1:p) {
  params[m.train + 1 + pp]  <- paste0("l", pp)
}
params[m.train + p + 2] <- "lp.f"
sims <- mcmc_array(ns, nchains = nchains, params)

for (ic in 1:nchains)
  sims[, ic, ] <- t(gp.class(dist = D, Y = Y, split_index = split_index
                             , a = a, b = b, ns = ns, monitor = FALSE))
pred.GP_lambda <- gp.class_pred(sims[, 1 ,], D, split_index, p
                                , burn = .2, thin = 10, avg = TRUE)
occ.GP_lambda <- gp.occ(sims[,1,], D, burn = .2, thin = 10)

# GP-RW ----
# GP with random walk kernel
X <- lapply(G, graph_from_adjacency_matrix, mode = "undirected", weighted = TRUE)
X <- lapply(X, function(x) x %>% set_vertex_attr("label", value = V(x)))
K <- CalculateKStepRandomWalkKernel(X, rep(1, 3))

ns <- 5000
nchains <- 1
m.train <- sum(split_index == "train")
params <- list(m.train + 1)
for (i in 1:m.train) {
  params[i] <- paste0("f", i)
}
params[m.train + 1] <- "lp.f"
sims <- mcmc_array(ns, nchains = nchains, params)

for (ic in 1:nchains)
  sims[, ic, ] <- t(gp_rw.class(K = K, Y = Y, split_index = split_index
                                , ns = ns, monitor = FALSE))
pred.GP_rw <- gp_rw.class_pred(sims[, 1 ,], K, split_index
                               , burn = .2, thin = 10, avg = TRUE)
occ.GP_rw <- gp_rw.occ(sims[,1,], K, burn = .2, thin = 10)

# Evaluate ----
ifelse(pred.GP_F > .5, 1, -1)
ifelse(pred.GP_lambda > .5, 1, -1)
ifelse(pred.GP_rw > .5, 1, -1)

# Visualize ----
# pi vs mu vs sigma vs heuristic
p_F <- plot.occ(occ.GP_F, y_lab = "GP-F", x_lab = FALSE)
p_lambda <- plot.occ(occ.GP_lambda, y_lab = "GP-λ", x_lab = FALSE)
p_rw <- plot.occ(occ.GP_rw, y_lab = "GP-RW")

grid.arrange(p_F, p_lambda, p_rw, nrow = 3)