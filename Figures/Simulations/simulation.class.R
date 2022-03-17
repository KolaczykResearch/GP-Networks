# Version : R 4.0.0
# Libraries ----
library(parallel)
library(ergm)
library(kernlab)
library(graphkernels)
library(class)
library(graphclass) # devtools::install_github("jesusdaniel/graphclass")

# Functions ----
source(".R/rcorrER.R")
source(".R/gp_functions.gibbs_slice.R")
source(".R/revisions/dist.R")
source(".R/gp_rw.gibbs.R")
source(".R/mase.R")

compare_classifiers <- function(m, n, mod) {
  
  # Setting ----
  if (mod == 1) { # small-world
    p0 <- .05
    p1 <- .07
    G0 <- replicate(m/2, as_adjacency_matrix(simplify(sample_smallworld(1, n, 5, p0))), simplify = FALSE)
    G1 <- replicate(m/2, as_adjacency_matrix(simplify(sample_smallworld(1, n, 5, p1))), simplify = FALSE)    
    
  } else if (mod == 2) { # SBM
    pm0 <- cbind(c(.05, .15)
                 , c(.15, .05))
    pm1 <- cbind(c(.1, .15)
                 , c(.15, .05))
    G0 <- replicate(m/2, as_adj(sample_sbm(n, pm0, c(n/2, n/2))))
    G1 <- replicate(m/2, as_adj(sample_sbm(n, pm1, c(n/2, n/2))))
    
  } else if (mod == 3) { # corrER
    p0 <- p1 <- .8
    r0 <- r1 <- .8
    # generate parent adjacency matrix
    A0 <- erdos.renyi.game(n, p0, type = "gnp")
    A1 <- erdos.renyi.game(n, p1, type = "gnp")
    # sample corrER networks
    G0 <- rcorrER(m/2, A0, r0)
    G1 <- rcorrER(m/2, A1, r1)
    
  } else if (mod == 4) { # preferential attachment
    p0 <- .6
    p1 <- 1.4
    G0 <- replicate(m/2, as_adj(sample_pa(n, power = p0, directed = FALSE)))
    G1 <- replicate(m/2, as_adj(sample_pa(n, power = p1, directed = FALSE)))
    
  } else if (mod == 5) { # ERGM
    tmp <- network(n, directed = FALSE, density = 0)
    
    beta0 <- -1
    beta1 <- 1
    G0 <- replicate(m/2, as.matrix.network.adjacency(simulate(tmp ~ edges + kstar(2) + triangle
                                                              , coef = c(0, -1, beta0)))
                    , simplify = FALSE)
    G1 <- replicate(m/2, as.matrix.network.adjacency(simulate(tmp ~ edges + kstar(2) + triangle
                                                              , coef = c(0, -1, beta1)))
                    , simplify = FALSE)
    
  }
  
  G <- c(G0, G1)
  
  Y <- c(rep(-1, m/2), rep(1, m/2))
  
  # Split data
  train_percent <- 6/8
  validate_percent <- 1 - 6/8
  split <- c(train = train_percent, validate = validate_percent)
  split_index <- sample(cut(seq(m)
                            , m * cumsum(c(0, split))
                            , labels = names(split)))
  
  # GP-F ----
  # GP with squared-exponential kernel and Frobenius distance
  D <- dist.frobenius(list(G))
  a <- c(10, 10)
  b <- c(50, 50)
  
  ns <- 1000
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
  
  # GP-Lambda ----
  # GP with squared-exponential kernel and spectral-based distance
  D <- dist.eigen(list(G), normalized = TRUE)
  a <- c(10, 10)
  b <- c(50, 50)
  
  ns <- 1000
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
  pred.GP_lambda <- gp.class_pred(sims[, 1 ,], D, split_index, p
                                  , burn = .2, thin = 10, avg = TRUE)
  
  # RW kernel ----
  X <- lapply(G, graph_from_adjacency_matrix, mode = "undirected", weighted = TRUE)
  if (mod != 5) {
    X <- lapply(X, function(x) x %>% set_vertex_attr("label", value = V(x)))    
  }
  X <- CalculateKStepRandomWalkKernel(X, rep(1, 3))
  
  # GP-RW ----
  # GP with random walk kernel
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
    sims[, ic, ] <- t(gp_rw.class(K = X, Y = Y, split_index = split_index
                                  , ns = ns, monitor = FALSE))
  pred.GP_rw <- gp_rw.class_pred(sims[, 1 ,], X, split_index
                                 , burn = .2, thin = 10, avg = TRUE)
  
  # SVM-RW ----
  # SVM with random walk kernel
  X_train <- X[split_index == "train", split_index == "train"]
  Y_train <- Y[split_index == "train"]
  
  cost.cv <- simplify2array(
    lapply(rev(2^(seq(-5, 15, 4)))
           , function(cost) {
             fit <- ksvm(X_train, Y_train, type = "C-svc", kernel = "matrix"
                         , C = cost, cross = 5)
             return(list(cost, fit@cross))
           })
  )
  
  cost.best <- cost.cv[, which.min(cost.cv[2, ])]
  f_svm <- ksvm(X_train, Y_train, type = "C-svc", kernel = "matrix"
                , C = cost.best[[1]]
                , prob.model = TRUE)
  
  X_validate <- as.kernelMatrix(X[split_index == "validate"
                                  , split_index == "train"][, SVindex(f_svm)
                                                            , drop = FALSE])
  Y_validate <- Y[split_index == "validate"]
  pred.SVM_RW <- predict(f_svm, X_validate)
  prob.SVM_RW <- try(predict(f_svm, X_validate, type = "prob")[, 1])
  if (class(prob.SVM_RW) == "try-error") { # line search fails
    prob.SVM_RW <- NA 
  }
  
  # graphclass ----
  X <- G
  X <- lapply(X, as.matrix)
  X <- t(sapply(X, function(x) x[upper.tri(x)]))
  gc <- graphclass(X = X[split_index == "train", ]
                   , Y = Y_train
                   , Xtest = X[split_index == "validate", ]
                   , Ytest = Y_validate)
  
  # MASE ----
  # multiple adjacency spectral embedding with kNN
  mase.embed <- try(mase(G[split_index == "train"], d = 2))
  if (class(mase.embed) == "try-error") { # eigen decomposition fails sometimes
    pred.mase <- prob.mase <- NA
  } else {
    mase.validate <- project_networks(G[split_index == "validate"], mase.embed$V)
    mase.Xtrain <- t(sapply(mase.embed$R, as.vector))
    mase.Xvalidate <- t(sapply(mase.validate, as.vector))
    pred.mase <- knn(train = mase.Xtrain, test = mase.Xvalidate
                     , cl = Y_train, k = 7, prob = TRUE)
    prob.mase <- attributes(pred.mase)$prob
  }
  
  # Evaluate ----
  m.validate <- sum(split_index == "validate")
  return(list(GP_F = sum(ifelse(pred.GP_F > .5, 1, -1) == Y_validate) / m.validate
              , GP_F.auc = roc(Y_validate, pred.GP_F)
              , GP_lambda = sum(ifelse(pred.GP_lambda > .5, 1, -1) == Y_validate) / m.validate
              , GP_lambda.auc = roc(Y_validate, pred.GP_lambda)
              , GP_RW = sum(ifelse(pred.GP_rw > .5, 1, -1) == Y_validate) / m.validate
              , GP_RW.auc = roc(Y_validate, pred.GP_rw)
              , SVM_RW = sum(pred.SVM_RW == Y_validate) / m.validate
              , SVM_RW.auc = roc(Y_validate, prob.SVM_RW)
              , GC = 1 - gc$test_error
              , GC.auc = roc(Y_validate, gc$Yfit_test)
              , MASE = sum(pred.mase == Y_validate) / m.validate
              , MASE.auc = roc(Y_validate, prob.mase)
              , m = m
              , n = n))
  
}

# Simulation ----
mod <- 1:5
m <- seq(40, 340, 60) # total sample size
n <- seq(20, 100, 20) # number of nodes
grid <- expand.grid(m, n, mod, stringsAsFactors = FALSE)
reps <- 100

results <- matrix(nrow = nrow(grid), ncol = 27)
results <- as.data.frame(results)
colnames(results) <- c("m", "n", "mod"
                       , "GP_F.acc", "GP_F.sd", "GP_F.auc", "GP_F.auc_sd"
                       , "GP_lambda.acc", "GP_lambda.sd", "GP_lambda.auc", "GP_lambda.auc_sd"
                       , "GP_RW.acc", "GP_RW.sd", "GP_RW.auc", "GP_RW.auc_sd"
                       , "SVM_RW.acc", "SVM_RW.sd", "SVM_RW.auc", "SVM_RW.auc_sd"
                       , "GC.acc", "GC.sd", "GC.auc", "GC.auc_sd"
                       , "MASE.acc", "MASE.sd", "MASE.auc", "MASE.auc_sd")

n.cores <- detectCores()
for (g in seq_len(nrow(grid))) {
  # Set simulation condition
  m <- grid[g, 1]
  n <- grid[g, 2]
  mod <- grid[g, 3]
  
  # Replicate condition
  start <- Sys.time()
  result <- simplify2array(mclapply(seq_len(reps), mc.cores = n.cores
                                    , FUN = function (r) {
                                      set.seed(r) # reproducibility in parallel
                                      unlist(compare_classifiers(m = m, n = n, mod = mod))
                                    }))
  message(g, " (", round(difftime(Sys.time(), start, units = "secs")), " secs) : "
          , " m = ", m
          , " , n = ", n
          , " , mod = ", mod)
  
  GP_F <- result["GP_F", ]
  GP_lambda <- result["GP_lambda", ]
  GP_RW <- result["GP_RW", ]
  SVM_RW <- result["SVM_RW", ]
  GC <- result["GC", ]
  MASE <- result["MASE", ]
  
  GP_F.auc <- result["GP_F.auc", ]
  GP_lambda.auc <- result["GP_lambda.auc", ]
  GP_RW.auc <- result["GP_RW.auc", ]
  SVM_RW.auc <- result["SVM_RW.auc", ]
  GC.auc <- result["GC.auc", ]
  MASE.auc <- result["MASE.auc", ]
  
  results[g, ] <- c(m, n, mod
                    , mean(GP_F), sd(GP_F), mean(GP_F.auc), sd(GP_F.auc)
                    , mean(GP_lambda), sd(GP_lambda), mean(GP_lambda.auc), sd(GP_lambda.auc)
                    , mean(GP_RW), sd(GP_RW), mean(GP_RW.auc), sd(GP_RW.auc)
                    , mean(SVM_RW), sd(SVM_RW), mean(SVM_RW.auc, na.rm = TRUE), sd(SVM_RW.auc, na.rm = TRUE)
                    , mean(GC), sd(GC), mean(GC.auc), sd(GC.auc)
                    , mean(MASE, na.rm = TRUE), sd(MASE, na.rm = TRUE), mean(MASE.auc, na.rm = TRUE), sd(MASE.auc, na.rm = TRUE))
}

# Plot ----
library(reshape2)
library(ggplot2)

# accuracy
df.acc <- melt(results[, c(1:3, seq(4, 24, 4))], id = c("m", "n", "mod"))
colnames(df.acc) <- c("m", "n", "Setting", "Method", "Accuracy")
df.acc$Method <- sub("\\..*", "", df.acc[, "Method"])
df.sd <- melt(results[, c(1:3, seq(5, 25, 4))], id = c("m", "n", "mod"))
colnames(df.sd) <- c("m", "n", "Setting", "Method", "sd")
df.sd$Method <- sub("\\..*", "", df.sd[, "Method"])
df <- merge(df.acc, df.sd)

df$n <- factor(df$n, levels = sort(unique(df$n)), labels = paste0("n = ", sort(unique(df$n))))
df$Method <- factor(df$Method, levels = c("GP_F", "GP_lambda", "GP_RW", "SVM_RW", "GC", "MASE"), labels = c("GP-F", "GP-λ", "GP-RW", "SVM-RW", "GC", "MASE"))
df$Setting <- factor(df$Setting, levels = 1:5, labels = c("small-world", "SBM", "corrER", "PA", "ERGM"))

df$Accuracy <- as.numeric(df$Accuracy)
df$sd <- as.numeric(df$sd)

ggplot(df, aes(x = m, y = Accuracy, group = Method, color = Method)) +
  # ylim(c(0, 1)) +
  geom_errorbar(data = df, aes(ymin = pmax(0, Accuracy - sd)
                               , ymax = pmin(1, Accuracy + sd)
                               , color = Method
                               , width = 60)
                , position = position_dodge(width = 10)) +
  geom_point(aes(shape = Method), position = position_dodge(width = 10)) +
  geom_line(position = position_dodge(width = 5)) +
  scale_x_continuous(breaks = unique(results$m)[c(1, 3, 5)]) +
  facet_grid(Setting ~ n) +
  theme_bw()

# AUC
df.auc <- melt(results[, c(1:3, seq(6, 26, 4))], id = c("m", "n", "mod"))
colnames(df.auc) <- c("m", "n", "Setting", "Method", "AUC")
df.auc$Method <- sub("\\..*", "", df.acc[, "Method"])
df.sd <- melt(results[, c(1:3, seq(7, 27, 4))], id = c("m", "n", "mod"))
colnames(df.sd) <- c("m", "n", "Setting", "Method", "sd")
df.sd$Method <- sub("\\..*", "", df.sd[, "Method"])
df <- merge(df.auc, df.sd)

df$n <- factor(df$n, levels = sort(unique(df$n)), labels = paste0("n = ", sort(unique(df$n))))
df$Method <- factor(df$Method, levels = c("GP_F", "GP_lambda", "GP_RW", "SVM_RW", "GC", "MASE"), labels = c("GP-F", "GP-λ", "GP-RW", "SVM-RW", "GC", "MASE"))
df$Setting <- factor(df$Setting, levels = 1:5, labels = c("small-world", "SBM", "corrER", "PA", "ERGM"))

df$AUC <- as.numeric(df$AUC)
df$sd <- as.numeric(df$sd)

ggplot(df, aes(x = m, y = AUC, group = Method, color = Method)) +
  # ylim(c(0, 1)) +
  geom_errorbar(data = df, aes(ymin = pmax(0, AUC - sd)
                               , ymax = pmin(1, AUC + sd)
                               , color = Method
                               , width = 60)
                , position = position_dodge(width = 10)) +
  geom_point(aes(shape = Method), position = position_dodge(width = 10)) +
  geom_line(position = position_dodge(width = 5)) +
  scale_x_continuous(breaks = unique(results$m)[c(1, 3, 5)]) +
  facet_grid(Setting ~ n) +
  theme_bw()
