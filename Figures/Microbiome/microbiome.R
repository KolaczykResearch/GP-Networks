# Version : R 4.0.5
# Libraries ----
library(ggplot2)
library(grid)
library(dplyr)
library(foreach)
library(doParallel)
library(graphkernels)

# Data ----
load(".Data/microbiome.Rdata")
m <- length(G)

# Classification ----
source("./R/gp_functions.gibbs_slice.R")
source("./R/revisions/dist.R")
source("./R/gp_rw.gibbs.R")

# GP-F ----
D_F <- dist.frobenius(list(G))
a <- c(0, 10)
b <- c(0, 50)

# Pre-allocate
ns <- 10000
nchains <- 1
m.train <- m - 1 # LOOCV
p <- dim(D_F)[3]
params <- list(m.train + p + 2)
for (i in seq_len(m.train)) {
  params[i] <- paste0("f", i)
}
params[m.train + 1] <- "sigma"
for (pp in seq_len(p)) {
  params[m.train + 1 + pp]  <- paste0("l", pp)
}
params[m.train + p + 2] <- "lp.f"

cores <- detectCores()
cl <- makeCluster(cores)
registerDoParallel(cl)
y.prob_F <- foreach(i = seq_len(m), .combine = rbind, .export = "slice"
                    , .packages = "invgamma"
) %dopar% {
  # train/test split
  split_index <- rep("train", m)
  split_index[i] <- "validate"
  
  set.seed(575)
  sims <- mcmc_array(ns, nchains = nchains, params)
  for (ic in seq_len(nchains))
    sims[, ic, ] <- t(gp.class(dist = D_F, Y = Y, split_index = split_index
                               , a = a, b = b, ns = ns, monitor = FALSE))
  
  # mixing
  # plot(sims[,1,][, m.train + 1], type = "l")    # sigma
  # plot(sims[,1,][, m.train + 2], type = "l")    # ell
  # plot(sims[,1,][, 1], type = "l")              # f
  
  # predictions
  return(gp.class_pred(sims[,1,], D_F, split_index, p
                       , burn = .2, thin = 10, avg = TRUE))
} # end for loop
stopCluster(cl)

roc(Y, y.prob_F, plot = TRUE)
y.pred_F <- ifelse(y.prob_F > .5, 1, ifelse(y.prob_F < .5, -1, NA))
table(Y, y.pred_F)
mcc(Y, y.pred_F)
f1(Y, y.pred_F)

# rerun with AVG = FALSE for prediction intervals
# y.pred_F <- ifelse(y.prob_F[, 1] > .5, 1, ifelse(y.prob_F[, 1] < .5, -1, NA))
# y.pred_F.sig <- ifelse(sapply(seq_len(m), function(i) .5 >= y.prob_F[i, 2]
#                               & .5 <= y.prob_F[i, 3]), 0, 1)
# table(Y, y.pred_F*y.pred_F.sig)

# GP-lambda ----
D_lambda <- dist.eigen(list(G), normalized = TRUE, signed = TRUE)

# Pre-allocate
ns <- 10000
nchains <- 1
m.train <- m - 1 # LOOCV
p <- dim(D_lambda)[3]
params <- list(m.train + p + 2)
for (i in seq_len(m.train)) {
  params[i] <- paste0("f", i)
}
params[m.train + 1] <- "sigma"
for (pp in seq_len(p)) {
  params[m.train + 1 + pp]  <- paste0("l", pp)
}
params[m.train + p + 2] <- "lp.f"

cores <- detectCores()
cl <- makeCluster(cores)
registerDoParallel(cl)
y.prob_lambda <- foreach(i = seq_len(m), .combine = rbind, .export = "slice"
                         , .packages = "invgamma"
) %dopar% {
  # train/test split
  split_index <- rep("train", m)
  split_index[i] <- "validate"
  
  set.seed(575)
  sims <- mcmc_array(ns, nchains = nchains, params)
  for (ic in seq_len(nchains))
    sims[, ic, ] <- t(gp.class(dist = D_lambda, Y = Y, split_index = split_index
                               , a = a, b = b, ns = ns, monitor = FALSE))
  
  # mixing
  # plot(sims[,1,][, m.train + 1], type = "l")    # sigma
  # plot(sims[,1,][, m.train + 2], type = "l")    # ell
  # plot(sims[,1,][, 1], type = "l")              # f
  
  # predictions
  return(gp.class_pred(sims[,1,], D_lambda, split_index, p
                       , burn = .2, thin = 10, avg = TRUE))
} # end for loop
stopCluster(cl)

roc(Y, y.prob_lambda, plot = TRUE)
y.pred_lambda <- ifelse(y.prob_lambda > .5, 1, ifelse(y.prob_lambda < .5, -1, NA))
table(Y, y.pred_lambda)
mcc(Y, y.pred_lambda)
f1(Y, y.pred_lambda)

# rerun with AVG = FALSE for prediction intervals
# y.pred_lambda <- ifelse(y.prob_lambda[, 1] > .5, 1, ifelse(y.prob_lambda[, 1] < .5, -1, NA))
# y.pred_lambda.sig <- ifelse(sapply(seq_len(m), function(i) .5 >= y.prob_lambda[i, 2]
#                               & .5 <= y.prob_lambda[i, 3]), 0, 1)
# table(Y, y.pred_lambda*y.pred_lambda.sig)

# GP-RW ----
G_bin <- lapply(G, function(A) ifelse(abs(A) > .1, 1, 0))
mean(sapply(G_bin, function(g) sum(g[lower.tri(g)])))

K <- lapply(G_bin, graph_from_adjacency_matrix, mode = "undirected", weighted = TRUE)
K <- lapply(K, function(x) x %>% set_vertex_attr("label", value = V(x)))
K <- CalculateKStepRandomWalkKernel(K, rep(1, 3))

# Pre-allocate
ns <- 50000
nchains <- 1
m.train <- m - 1 # LOOCV
params <- list(m.train + 1)
for (i in seq_len(m.train)) {
  params[i] <- paste0("f", i)
}
params[m.train + 1] <- "lp.f"

cores <- detectCores()
cl <- makeCluster(cores)
registerDoParallel(cl)
y.prob_rw <- foreach(i = seq_len(m), .combine = rbind, .export = "slice") %dopar% {
  # train/test split
  split_index <- rep("train", m)
  split_index[i] <- "validate"
  
  set.seed(575)
  sims <- mcmc_array(ns, nchains = nchains, params)
  
  for (ic in seq_len(nchains))
    sims[, ic, ] <- t(gp_rw.class(K, Y = Y, split_index = split_index
                                  , ns = ns))
  
  # predictions
  return(gp_rw.class_pred(sims[, 1 ,], K, split_index
                          , burn = .4, thin = 10, avg = TRUE))
} # end for loop
stopCluster(cl)

roc(Y, y.prob_rw, plot = TRUE)
y.pred_rw <- ifelse(y.prob_rw > .5, 1, ifelse(y.prob_rw < .5, -1, NA))
table(Y, y.pred_rw)
mcc(Y, y.pred_rw)
f1(Y, y.pred_rw)

# rerun with AVG = FALSE for prediction intervals
# y.pred_rw <- ifelse(y.prob_rw[, 1] > .5, 1, ifelse(y.prob_rw[, 1] < .5, -1, NA))
# y.pred_rw.sig <- ifelse(sapply(seq_len(m), function(i) .5 >= y.prob_rw[i, 2]
#                               & .5 <= y.prob_rw[i, 3]), 0, 1)
# table(Y, y.pred_rw*y.pred_rw.sig)

# OCC ----
set.seed(575)
source("./R/gp_functions.occ.R")

split_index <- ifelse(Y == -1, "train", "validate")
ind.train <- which(split_index == "train")
ind.validate <- which(split_index == "validate")
split_index[sample(ind.validate, 3)] <- "train"
split_index[sample(ind.train, 5)] <- "validate"

# GP-F ----
a <- c(0, 1)
b <- c(0, 1)

ns <- 10000
nchains <- 1
m.train <- sum(split_index == "train")
p <- dim(D_F)[3]
params <- list(m.train + p + 2)
for (i in seq_len(m.train)) {
  params[i] <- paste0("f", i)
}
params[m.train + 1] <- "sigma"
for (pp in seq_len(p)) {
  params[m.train + 1 + pp]  <- paste0("l", pp)
}
params[m.train + p + 2] <- "lp.f"

sims <- mcmc_array(ns, nchains = nchains, params)
for (ic in seq_len(nchains))
  sims[, ic, ] <- t(gp.class(dist = D_F, Y = Y, split_index = split_index
                             , a = a, b = b, ns = ns, monitor = FALSE))
occ.GP_F <- gp.occ(sims[,1,], D_F, burn = .2, thin = 10)
sapply(occ.GP_F, function(score) roc(labels = Y[split_index == "validate"]
                                     , scores = score))

# GP-lambda ----
ns <- 10000
nchains <- 1
m.train <- sum(split_index == "train")
p <- dim(D_lambda)[3]
params <- list(m.train + p + 2)
for (i in seq_len(m.train)) {
  params[i] <- paste0("f", i)
}
params[m.train + 1] <- "sigma"
for (pp in seq_len(p)) {
  params[m.train + 1 + pp]  <- paste0("l", pp)
}
params[m.train + p + 2] <- "lp.f"

sims <- mcmc_array(ns, nchains = nchains, params)
for (ic in seq_len(nchains))
  sims[, ic, ] <- t(gp.class(dist = D_lambda, Y = Y, split_index = split_index
                             , a = a, b = b, ns = ns, monitor = FALSE))
occ.GP_lambda <- gp.occ(sims[,1,], D_lambda, burn = .2, thin = 10)
sapply(occ.GP_lambda, function(score) roc(labels = Y[split_index == "validate"]
                                          , scores = score))

save.image(file = "./Data/occ_lambda.RData")

# GP-RW ----
ns <- 50000
nchains <- 1
m.train <- sum(split_index == "train")
params <- list(m.train + 1)
for (i in seq_len(m.train)) {
  params[i] <- paste0("f", i)
}
params[m.train + 1] <- "lp.f"
sims <- mcmc_array(ns, nchains = nchains, params)

for (ic in seq_len(nchains))
  sims[, ic, ] <- t(gp_rw.class(K = K, Y = Y, split_index = split_index
                                , ns = ns, monitor = FALSE))
occ.GP_rw <- gp_rw.occ(sims[,1,], K, burn = .2, thin = 10)
sapply(occ.GP_rw, function(score) roc(labels = Y[split_index == "validate"]
                                      , scores = score))

p_F <- plot.occ(occ.GP_F, y_lab = "GP-F", x_lab = FALSE)
p_lambda <- plot.occ(occ.GP_lambda, y_lab = "GP-λ", x_lab = FALSE)
p_rw <- plot.occ(occ.GP_rw, y_lab = "GP-RW")

grid.arrange(p_F, p_lambda, p_rw, nrow = 3)

# Survival ----
source("./R/gp_functions.survival.R")

Y <- GDDel
D_lambda <- dist.eigen(list(G), cbind(X, Y))
a <- c(1, 0, 10, 10, 10)
b <- c(0, 0, 50, 5, 5)

ns <- 1000

set.seed(575)
system.time(surv.lambda <- gp.survival(D_lambda, Y, a, b, N = ns))

# plot(surv.lambda$Omega, typ = "l")            # Omega
# plot(surv.lambda$hyperparams[1, ], typ = "l") # sigma
# plot(surv.lambda$hyperparams[2, ], typ = "l") # ell.net
# plot(surv.lambda$hyperparams[3, ], typ = "l") # ell.X1
# plot(surv.lambda$hyperparams[4, ], typ = "l") # ell.time
# plot(sapply(surv.lambda$l, function(l) l[1]), typ = "l") # l(1)
# plot(sapply(1:dim(surv.lambda$surv)[3]
#             , function(q) surv.lambda$surv[1, 1, q]), typ = "l") # l(grid(1), 1)

p1 <- plot.surv(surv.lambda, Y, groups = X, km = T)

p1

# add to plot.surv

# switch names in p1
# geom_line(alpha = .2)
# scale_colour_manual(name = "")

# p1 <- p1 +
#   scale_x_continuous(breaks = c(0, 10, 20, 30, 40) * 7 / 100
#                      , labels = c(0, 10, 20, 30, 40)) +
#   xlab("Gestational Weeks")