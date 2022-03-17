# Version : R 4.0.5
# Libraries ----
library(parallel)
library(MASS)

# Functions ----
source(".R/dist.R")
source(".R/gp_functions.gibbs_slice.R")

compare_classifiers <- function(m, n, n_k) {
  
  # Sample data ----
  # sample 2 covariance matrices
  A0 <- matrix(runif(n^2) - .5, ncol = n) # matrix(runif(n^2)*2-1, ncol = n)
  Sigma0 <- crossprod(A0)
  
  A1 <- matrix(runif(n^2) - .5, ncol = n)
  Sigma1 <- crossprod(A1)
  
  # sample data from N(0, Sigma)
  X0 <- replicate(m/2, mvrnorm(n = n_k, mu = rep(0, n), Sigma = Sigma0)
                  , simplify = FALSE)
  X1 <- replicate(m/2, mvrnorm(n = n_k, mu = rep(0, n), Sigma = Sigma1)
                  , simplify = FALSE)
  X <- c(X0, X1)
  
  G0 <- lapply(X0, cor)
  G0 <- lapply(G0, function(A) {diag(A) <- 0; A})
  
  G1 <- lapply(X1, cor)
  G1 <- lapply(G1, function(A) {diag(A) <- 0; A})
  
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
  D <- dist.eigen(list(G), normalized = TRUE, signed = TRUE)
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
  
  # GP-raw ----
  # GP with squared-exponential kernel on non-network data
  D <- as.matrix(dist(t(sapply(X, as.vector)))^2) / (n*(n-1))
  D <- array(D, dim = c(m, m, 1))
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
  pred.GP_raw <- gp.class_pred(sims[, 1 ,], D, split_index, p
                               , burn = .2, thin = 10, avg = TRUE)
  
  # Evaluate ----
  m.validate <- sum(split_index == "validate")
  Y_validate <- Y[split_index == "validate"]
  return(list(GP_F = sum(ifelse(pred.GP_F > .5, 1, -1) == Y_validate) / m.validate
              , GP_F.auc = roc(Y_validate, pred.GP_F)
              , GP_lambda = sum(ifelse(pred.GP_lambda > .5, 1, -1) == Y_validate) / m.validate
              , GP_lambda.auc = roc(Y_validate, pred.GP_lambda)
              , GP_raw = sum(ifelse(pred.GP_raw > .5, 1, -1) == Y_validate) / m.validate
              , GP_raw.auc = roc(Y_validate, pred.GP_raw)
              , m = m
              , n = n))
  
}

# Simulation ----
m <- seq(20, 50, 10)  # total sample size (subjects)
n <- seq(20, 180, 40) # number of nodes (taxa)
n_k <- seq(5, 25, 10) # number of samples
reps <- 100

grid <- expand.grid(m, n, n_k, stringsAsFactors = FALSE)

results <- matrix(nrow = nrow(grid), ncol = 15)
results <- as.data.frame(results)
colnames(results) <- c("m", "n", "n_k"
                       , "GP_F.acc", "GP_F.sd", "GP_F.auc", "GP_F.auc_sd"
                       , "GP_lambda.acc", "GP_lambda.sd", "GP_lambda.auc", "GP_lambda.auc_sd"
                       , "GP_raw.acc", "GP_raw.sd", "GP_raw.auc", "GP_raw.auc_sd")

n.cores <- detectCores()
for (g in seq_len(nrow(grid))) {
  # Set simulation condition
  m <- grid[g, 1]
  n <- grid[g, 2]
  n_k <- grid[g, 3]
  
  # Replicate condition
  start <- Sys.time()
  result <- simplify2array(mclapply(seq_len(reps), mc.cores = n.cores
                                    , FUN = function (r) {
                                      set.seed(r) # reproducibility in parallel
                                      unlist(compare_classifiers(m = m, n = n, n_k = n_k))
                                    }))
  message(g, " (", round(difftime(Sys.time(), start, units = "secs")), " secs) : "
          , " m = ", m
          , " , n = ", n
          , " , n_k = ", n_k)
  
  GP_F <- result["GP_F", ]
  GP_lambda <- result["GP_lambda", ]
  GP_raw <- result["GP_raw", ]
  
  GP_F.auc <- result["GP_F.auc", ]
  GP_lambda.auc <- result["GP_lambda.auc", ]
  GP_raw.auc <- result["GP_raw.auc", ]
  
  results[g, ] <- c(m, n, n_k
                    , mean(GP_F), sd(GP_F), mean(GP_F.auc), sd(GP_F.auc)
                    , mean(GP_lambda), sd(GP_lambda), mean(GP_lambda.auc), sd(GP_lambda.auc)
                    , mean(GP_raw), sd(GP_raw), mean(GP_raw.auc), sd(GP_raw.auc))
}

# Plot ----
library(reshape2)
library(ggplot2)
library(scales)

# accuracy
df.acc <- melt(results[, c(1:3, seq(4, 12, 4))], id = c("m", "n", "n_k"))
colnames(df.acc) <- c("m", "n", "n_k", "Method", "Accuracy")
df.acc$Method <- sub("\\..*", "", df.acc[, "Method"])
df.sd <- melt(results[, c(1:3, seq(5, 13, 4))], id = c("m", "n", "n_k"))
colnames(df.sd) <- c("m", "n", "n_k", "Method", "sd")
df.sd$Method <- sub("\\..*", "", df.sd[, "Method"])
df <- merge(df.acc, df.sd)

df$n <- factor(df$n, levels = sort(unique(df$n)), labels = c(expression(n[]*" = 20"), expression(n[]*" = 60"), expression(n[]*" = 100"), expression(n[]*" = 140"), expression(n[]*" = 180")))
df$Method <- factor(df$Method, levels = c("GP_F", "GP_lambda", "GP_raw"), labels = c("GP-F", "GP-λ", "GP-raw"))
df$n_k <- factor(df$n_k, levels = sort(unique(df$n_k)), labels = c(expression(N[k]*" = 5"), expression(N[k]*" = 15"), expression(N[k]*" = 25")))

df$Accuracy <- as.numeric(df$Accuracy)
df$sd <- as.numeric(df$sd)

# for consistency with other figure
cols <- hue_pal()(6)[c(1, 2)]
cols <- c(cols, "black")

ggplot(df, aes(x = m, y = Accuracy, group = Method, color = Method)) +
  # ylim(c(0, 1)) +
  geom_errorbar(data = df, aes(ymin = pmax(0, Accuracy - sd)
                               , ymax = pmin(1, Accuracy + sd)
                               , color = Method
                               , width = 10)
                , position = position_dodge(width = 5)) +
  geom_point(aes(shape = Method), position = position_dodge(width = 5)) +
  geom_line(position = position_dodge(width = 5)) +
  scale_x_continuous(breaks = unique(results$m)) +
  facet_grid(n_k ~ n, labeller = label_parsed) +
  theme_bw() +
  scale_color_manual(values = cols)