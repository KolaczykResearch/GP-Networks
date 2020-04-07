# Libraries ----
library(foreach)
library(doParallel)

# Functions ----
source("/R/gp_functions.gibbs_slice.R")

# Data ----
# from : https://github.com/jesusdaniel/graphclass/tree/master/data
# processed using : get_matrix()
load("/Data/cobre.RData")
N <- length(Y)

# 10-fold CV ----
N.fold <- 10
cv.folds <- (1:N %% N.fold) + 1

# Inputs ----
D <- kern.dis((list(G)))

a <- c(0, 10) # flat prior on sigma
b <- c(0, 5)

cores <- detectCores()
cl <- makeCluster(cores[1] - 1)
registerDoParallel(cl)
cv <- foreach(fold = 1:N.fold, .combine = rbind, .multicombine = TRUE
              , .packages = "invgamma", .export = "slice"
) %dopar% {
  
  set.seed(575) # inside parallelization for reproducibility
  
  # train/test split
  split_index <- rep("train", N)
  split_index[which(cv.folds == fold)] <- "validate"
  
  # MCMC
  ns <- 50000
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
                               , a = a, b = b, ns = ns))
  
  # mixing
  # plot(sims[, 1 ,][, N.train + 1], type = "l") # sigma
  # plot(sims[, 1 ,][, N.train + 2], type = "l") # ell
  # plot(sims[, 1 ,][, 1], type = "l")           # f
  
  # predictions
  y.pred <- gp.class_pred(sims[, 1 ,], D, split_index, p
                          , burn = .4, thin = 10, avg = TRUE)
  
  c(sum(diag(prop.table(table(Y[split_index == "validate"]            # acc
                              , ifelse(y.pred > .5, 1
                                       , ifelse(y.pred < .5, -1, NA))))))
    , mcc(Y[split_index == "validate"]                                # mcc
          , ifelse(y.pred > .5, 1, ifelse(y.pred < .5, -1, NA)))
    , roc(Y[split_index == "validate"], y.pred))                      # auc
  
} # end for loop
stopCluster(cl)

# Results ----
cv <- as.data.frame(cv)
names(cv) <- c("acc", "mcc", "auc")
(acc.mean <- mean(cv$acc)) # COBRE : Acc = 91.9 (9.6)  ; AUC = 97.5 (4.7)
acc.sd <- sd(cv$acc)
mcc.mean <- mean(cv$mcc)
(auc.mean <- mean(cv$auc))