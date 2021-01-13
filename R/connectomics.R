# Libraries ----
library(foreach)
library(doParallel)
library(graphkernels)

# Functions ----
source("/R/gp_functions.gibbs_slice.R")
source("/R/dist.R")
source("/R/gp_rw.gibbs.R")

# Data ----
# from : https://github.com/jesusdaniel/graphclass/tree/master/data
# processed using : get_matrix()
load("/Data/cobre.RData")
m <- length(Y)

# 10-fold CV ----
N.fold <- 10
cv.folds <- (1:m %% N.fold) + 1

# GP-F ----
D <- dist.frobenius(list(G))
a <- c(0, 10)
b <- c(0, 5)

cores <- detectCores()
cl <- makeCluster(cores)
registerDoParallel(cl)
cv.GP_F <- foreach(fold = 1:N.fold, .combine = rbind, .multicombine = TRUE
                   , .packages = "invgamma", .export = "slice"
) %dopar% {
  
  set.seed(575) # inside parallelization for reproducibility
  
  # train/test split
  split_index <- rep("train", m)
  split_index[which(cv.folds == fold)] <- "validate"
  
  # MCMC
  ns <- 50000
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
                               , a = a, b = b, ns = ns))
  
  # mixing
  # plot(sims[, 1 ,][, m.train + 1], type = "l") # sigma
  # plot(sims[, 1 ,][, m.train + 2], type = "l") # ell
  # plot(sims[, 1 ,][, 1], type = "l")           # f
  
  # predictions
  pred.GP_F <- gp.class_pred(sims[, 1 ,], D, split_index, p
                             , burn = .4, thin = 10, avg = TRUE)
  
  m.validate <- sum(split_index == "validate")
  c(sum(ifelse(pred.GP_F > .5, 1, -1) == Y[split_index == "validate"]) / m.validate
    , mcc(Y[split_index == "validate"]
          , ifelse(pred.GP_F > .5, 1, -1))
    , roc(Y[split_index == "validate"], pred.GP_F))
  
} # end for loop
stopCluster(cl)

cv.GP_F <- as.data.frame(cv.GP_F)
names(cv.GP_F) <- c("acc", "mcc", "auc")
(acc.GP_F <- mean(cv.GP_F$acc)); acc_sd.GP_F <- sd(cv.GP_F$acc)
mcc.GP_F <- mean(cv.GP_F$mcc)
(auc.GP_F <- mean(cv.GP_F$auc)); auc_sd.GP_F <- sd(cv.GP_F$auc)

# GP-lambda ----
D <- dist.eigen(list(G), normalized = TRUE, signed = TRUE)
a <- c(0, 10)
b <- c(0, 5)

cores <- detectCores()
cl <- makeCluster(cores)
registerDoParallel(cl)
cv.GP_lambda <- foreach(fold = 1:N.fold, .combine = rbind, .multicombine = TRUE
                        , .packages = "invgamma", .export = "slice"
) %dopar% {
  
  set.seed(575) # inside parallelization for reproducibility
  
  # train/test split
  split_index <- rep("train", m)
  split_index[which(cv.folds == fold)] <- "validate"
  
  # MCMC
  ns <- 50000
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
                               , a = a, b = b, ns = ns))
  
  # mixing
  # plot(sims[, 1 ,][, m.train + 1], type = "l") # sigma
  # plot(sims[, 1 ,][, m.train + 2], type = "l") # ell
  # plot(sims[, 1 ,][, 1], type = "l")           # f
  
  # predictions
  pred.GP_lambda <- gp.class_pred(sims[, 1 ,], D, split_index, p
                                  , burn = .4, thin = 10, avg = TRUE)
  
  m.validate <- sum(split_index == "validate")
  c(sum(ifelse(pred.GP_lambda > .5, 1, -1) == Y[split_index == "validate"]) / m.validate
    , mcc(Y[split_index == "validate"]
          , ifelse(pred.GP_lambda > .5, 1, -1))
    , roc(Y[split_index == "validate"], pred.GP_lambda))
  
} # end for loop
stopCluster(cl)

cv.GP_lambda <- as.data.frame(cv.GP_lambda)
names(cv.GP_lambda) <- c("acc", "mcc", "auc")
(acc.GP_lambda <- mean(cv.GP_lambda$acc)); acc_sd.GP_lambda <- sd(cv.GP_lambda$acc)
mcc.GP_lambda <- mean(cv.GP_lambda$mcc)
(auc.GP_lambda <- mean(cv.GP_lambda$auc)); auc_sd.GP_lambda <- sd(cv.GP_lambda$auc)

# GP-RW ----
# Binarize :  ~5.3% of edges (mean 1833) because 5.4% (1886) reported as best
G_bin <- lapply(G, function(A) ifelse(abs(A) > .45, 1, 0))
mean(sapply(G_bin, function(g) sum(g[lower.tri(g)])))

K <- lapply(G_bin, graph_from_adjacency_matrix, mode = "undirected", weighted = TRUE)
K <- lapply(K, function(x) x %>% set_vertex_attr("label", value = V(x)))
K <- CalculateKStepRandomWalkKernel(K, rep(1, 2))

cores <- detectCores()
cl <- makeCluster(cores)
registerDoParallel(cl)
cv.GP_rw <- foreach(fold = 1:N.fold, .combine = rbind, .multicombine = TRUE
) %dopar% {
  
  set.seed(575) # inside parallelization for reproducibility
  
  # train/test split
  split_index <- rep("train", m)
  split_index[which(cv.folds == fold)] <- "validate"
  
  # MCMC
  ns <- 200000 # GP-RW is very fast because K is fixed
  nchains <- 1
  m.train <- sum(split_index == "train")
  params <- list(N.train + 1)
  for (i in 1:m.train) {
    params[i] <- paste0("f", i)
  }
  params[m.train + 1] <- "lp.f"
  sims <- mcmc_array(ns, nchains = nchains, params)
  
  for (ic in 1:nchains)
    sims[, ic, ] <- t(gp_rw.class(K, Y = Y, split_index = split_index
                                  , ns = ns))
  
  # predictions
  pred.GP_rw <- gp_rw.class_pred(sims[, 1 ,], K, split_index
                                 , burn = .4, thin = 10, avg = TRUE)
  
  m.validate <- sum(split_index == "validate")
  c(sum(ifelse(pred.GP_rw > .5, 1, -1) == Y[split_index == "validate"]) / m.validate
    , mcc(Y[split_index == "validate"]
          , ifelse(pred.GP_rw > .5, 1, -1))
    , roc(Y[split_index == "validate"], pred.GP_rw))
  
} # end for loop
stopCluster(cl)

cv.GP_rw <- as.data.frame(cv.GP_rw)
names(cv.GP_rw) <- c("acc", "mcc", "auc")
(acc.GP_rw <- mean(cv.GP_rw$acc)); acc_sd.GP_rw <- sd(cv.GP_rw$acc)
mcc.GP_rw <- mean(cv.GP_rw$mcc)
(auc.GP_rw <- mean(cv.GP_rw$auc)); auc_sd.GP_rw <- sd(cv.GP_rw$auc)

# Summary ----

# COBRE     :  Acc (sd)   ;  AUC (sd)
#
# GP-F      : .919 (.096)    .975 (.047)
# GP-lambda : .563 (.173)    .626 (.209)
# GP-RW     : .75  (.063)    .902 (.094)