library(parallel)
library(kernlab)

source("/R/gp_functions.gibbs_slice.R")

N <- c(20, 60, 100)
n <- c(10, 50, 100)
mod <- c("ERvsER", "SBMvsER", "SBMvsSBM")
dens <- c("low", "med", "high")
grid <- expand.grid(N, n, mod, dens, stringsAsFactors = FALSE)
rm(N, n, mod, dens)

sim <- function(N, n, mod, dens, ns = 2500) {
  
  if (mod == "ERvsER") {
    if (dens == "low") {
      p0 <- .1; p1 <- .2
    } else if (dens == "med") {
      p0 <- .45; p1 <- .55
    } else if (dens == "high") {
      p0 <- .85; p1 <- .95
    }
    G0 <- replicate(N/2, as_adj(erdos.renyi.game(n, p0, type = "gnp")), simplify = FALSE)
    G1 <- replicate(N/2, as_adj(erdos.renyi.game(n, p1, type = "gnp")), simplify = FALSE)
  } else if (mod == "SBMvsER") {
    if (dens == "low") {
      pm0 <- cbind(c(.15, .05)
                   , c(.05, .15))
      p <- .1
    } else if (dens == "med") {
      pm0 <- cbind(c(.55, .45)
                   , c(.45, .55))
      p <- .5
    } else if (dens == "high") {
      pm0 <- cbind(c(.95, .85)
                   , c(.85, .95))
      p <- .9
    }
    G0 <- replicate(N/2, as_adj(sample_sbm(n, pm0, c(n/2, n/2))))
    G1 <- replicate(N/2, as_adj(erdos.renyi.game(n, p, type = "gnp")), simplify = FALSE)
  } else if (mod == "SBMvsSBM") {
    if (dens == "low") {
      pm0 <- cbind(c(.15, .05)
                   , c(.05, .15))
      pm1 <- cbind(c(.25, .15)
                   , c(.15, .25))
    } else if (dens == "med") {
      pm0 <- cbind(c(.6, .5)
                   , c(.5, .6))
      pm1 <- cbind(c(.5, .4)
                   , c(.4, .5))
    } else if (dens == "high") {
      pm0 <- cbind(c(.9, .8)
                   , c(.8, .9))
      pm1 <- cbind(c(.8, .7)
                   , c(.7, .8))
    }
    G0 <- replicate(N/2, as_adj(sample_sbm(n, pm0, c(n/2, n/2))))
    G1 <- replicate(N/2, as_adj(sample_sbm(n, pm1, c(n/2, n/2))))
  }
  
  G <- c(G0, G1)                        
  Y <- c(rep(-1, N/2), rep(1, N/2))
  
  D <- kern.dis(list(G))
  a <- c(10, 10)
  b <- c(50, 50)
  
  # Split data
  train_percent <- 6/8
  validate_percent <- 1 - 6/8
  split <- c(train = train_percent, validate = validate_percent)
  split_index <- sample(cut(seq(N)
                            , N * cumsum(c(0, split))
                            , labels = names(split)))
  
  # GPC ----
  gpc.time <- Sys.time()
  ns <- ns
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
  
  # predictions
  y.pred <- gp.class_pred(sims[, 1 ,], D, split_index, p
                          , burn = .2, thin = 10, avg = TRUE)
  gpc.time <- Sys.time() - gpc.time
  
  # SVM ----
  svm.time <- Sys.time()
  svm.grid <- expand.grid(qinvgamma(seq(.05, .95, .1), a[1], b[1])    # sigma
                          , qinvgamma(seq(.05, .95, .1), a[2], b[2])) # ell
  
  K.cv <- simplify2array(
    lapply(1:nrow(svm.grid)
           , function(theta) {
             
             sigma <- svm.grid[theta, 1]; ell <- svm.grid[theta, 2]
             X <- as.kernelMatrix(kern.fun(D, c(sigma, ell)))
             
             X_train <- X[split_index == "train", split_index == "train"]
             Y_train <- Y[split_index == "train"]
             
             fit <- ksvm(X_train, Y_train, type = "C-svc", kernel = "matrix"
                         , C = 1, cross = 5)
             
             return(list(c(sigma, ell)
                         , X
                         , fit@cross))
           })
  )
  
  K.best <- K.cv[, which.min(K.cv[3, ])]
  
  X <- as.kernelMatrix(kern.fun(D, unlist(K.best[1])))
  X_train <- X[split_index == "train", split_index == "train"]
  Y_train <- Y[split_index == "train"]
  
  f_svm <- ksvm(X_train, Y_train, type = "C-svc", kernel = "matrix"
                , C = 1)
  
  X_validate <- as.kernelMatrix(X[split_index == "validate"
                                  , split_index == "train"][, SVindex(f_svm)
                                                            , drop = FALSE])
  Y_validate <- Y[split_index == "validate"]
  
  pred_validate <- predict(f_svm, X_validate)
  svm.time <- Sys.time() - svm.time
  
  return(list(GPC = table(Y[split_index == "validate"]
                          , ifelse(y.pred > .5, 1, ifelse(y.pred < .5, -1, NA)))
              , SVM = table(Y[split_index == "validate"], pred_validate)
              , GPC.time = gpc.time
              , SVM.time = svm.time
              , N = N
              , n = n
              , mod = mod
              , dens = dens))
}

reps <- 100
results <- matrix(nrow = nrow(grid), ncol = 10)
results <- as.data.frame(results)
colnames(results) <- c("N", "n", "mod", "dens"
                       , "GPC.acc", "GPC.sd", "GPC.time"
                       , "SVM.acc", "SVM.sd", "SVM.time")

for (g in 1:nrow(grid)) {
  print(g)
  
  # Set simulation condition
  N <- grid[g, 1]
  n <- grid[g, 2]
  mod <- grid[g, 3]
  dens <- grid[g, 4]
  
  # Replicate condition
  error <- simplify2array(mclapply(1:reps, mc.cores = 28, FUN = function (r) {
    set.seed(r) # reproducibility in parallel
    sim(N = N, n = n, mod = mod, dens = dens)
  }))
  
  gpc <- apply(error, 2, function(rep) sum(diag(prop.table(rep$GPC))))
  svm <- apply(error, 2, function(rep) sum(diag(prop.table(rep$SVM))))
  
  results[g, ] <- c(N, n, mod, dens
                    , mean(gpc), sd(gpc), mean(unlist(error["GPC.time", ]))
                    , mean(svm), sd(svm), mean(unlist(error["SVM.time", ])))
}

# Reshape ----

library(reshape2)
library(ggplot2)
library(gridExtra)

df.acc <- melt(results[, c(1:5, 8)], id = c("N", "n", "mod", "dens"))
colnames(df.acc) <- c("N", "n", "Model", "Density", "Classifier", "Accuracy")
df.acc$Classifier <- ifelse(df.acc$Classifier == "GPC.acc", "GPC", "SVM")
df.sd <- melt(results[, c(1:4, 6, 9)], id = c("N", "n", "mod", "dens"))
colnames(df.sd) <- c("N", "n", "Model", "Density", "Classifier", "sd")
df.sd$Classifier <- ifelse(df.sd$Classifier == "GPC.sd", "GPC", "SVM")
df <- merge(df.acc, df.sd)

df$N <- factor(df$N, levels = c(20, 60, 100))
df$n <- factor(df$n, levels = c(10, 50, 100))
levels(df$n) <- c("n = 10", "n = 50", "n = 100")
df$Density <- factor(df$Density, levels = c("high", "med", "low"))
levels(df$Density) <- c("High", "Medium", "Low")
df$Classifier <- factor(df$Classifier, levels = c("SVM", "GPC"))

df$Accuracy <- as.numeric(df$Accuracy)
df$sd <- as.numeric(df$sd)

# Figure 4 ----
ggplot(df[df$Model == "SBMvsER", ], aes(x = N, y = Accuracy, group = Classifier)) +
  ylim(c(0, 1)) +
  geom_errorbar(data = df[df$Model == "SBMvsER", ], aes(ymin = Accuracy - sd
                                                        , ymax = pmin(1, Accuracy + sd)
                                                        , color = Classifier
                                                        , width = .5)
                , position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5), color = "grey50") +
  facet_grid(Density ~ n) +
  theme_bw()

# Figure 10 ----
p1 <- ggplot(df[df$Model == "ERvsER", ], aes(x = N, y = Accuracy, group = Classifier)) +
  ylim(c(0, 1)) +
  geom_errorbar(data = df[df$Model == "ERvsER", ], aes(ymin = Accuracy - sd
                                                       , ymax = pmin(1, Accuracy + sd)
                                                       , color = Classifier
                                                       , width = .5)
                , position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5), color = "gray50") +
  facet_grid(Density ~ n) +
  theme_bw()

p2 <- ggplot(df[df$Model == "SBMvsSBM", ], aes(x = N, y = Accuracy, group = Classifier)) +
  ylim(c(0, 1)) +
  geom_errorbar(data = df[df$Model == "SBMvsSBM", ], aes(ymin = Accuracy - sd
                                                         , ymax = pmin(1, Accuracy + sd)
                                                         , color = Classifier
                                                         , width = .5)
                , position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5), color = "gray50") +
  facet_grid(Density ~ n) +
  theme_bw()

grid.arrange(p1, p2, ncol = 1)