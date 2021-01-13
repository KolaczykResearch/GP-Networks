# Libraries ----
library(igraph)
library(invgamma)

# Functions ----
kern.fun <- function(D, theta) {
  # computes squared-exponential kernel given distance matrix (array) D
  # NB: does not square distance
  
  N <- nrow(D)
  M <- ncol(D)
  p <- dim(D)[3]
  K0 <- matrix(0, N, M)
  sigma <- theta[1]
  
  for (pp in 1:p) {
    ell <- exp(theta[pp + 1])
    
    K0 <- K0 + exp(-ell * D[, , pp])
  }
  
  K <- sigma * K0  # sigma is overall signal variance
  
  return(K)
}

log1pe <- function (x) { # vectorized version: `x` can be a vector 
  l <- ifelse(x > 0, x, 0) # shift 
  x <- ifelse(x > 0, -x, x) # range reduction: `x = -abs(x)` 
  ifelse(x < log(.Machine$double.eps), l, l + log(1 + exp(x))) }

ilogit <- function (x) {
  if (x >= 0) {
    1 / (1 + exp(-x))  
  } else {
    z <- exp(x)
    z / (1 + z) 
  }
}

rmvnorm <- function(n = 1, mu, C){
  p <- length(mu)
  
  if (p == 1) {
    X <- rnorm(1, mu, C)
  } else{
    Z <- matrix(rnorm(p*n), p, n)
    
    X <- crossprod(C, Z)
    X <- sweep(X, 1, mu, FUN = `+`)
  }
  
  return(X)
}

chol.fun <- function(x, jit = 1e-6) {
  
  n <- nrow(x)
  L <- chol(x + diag(n) * jit^2)
  
  return(L)
}

slice <- function(theta, f, alpha, Y, D, a, b, sigma = 10) {
  
  log.lik <- function(u, f, C, g, theta, a, b, Y) {
    log(u) +                                                 # u
      -sum(log1pe(-Y*f)) +                                   # L(f)
      -.5 * (crossprod(backsolve(C, g, transpose = TRUE))) + # N(g; 0, Sigma + S)
      sum(-a * theta - b / exp(theta))                       # p(theta)
  }
  
  n <- length(f)
  p <- length(theta)
  Sigma <- kern.fun(D, c(alpha, theta))
  S <- diag(n) * alpha # auxillary noise
  S.inv <- diag(n) / alpha
  
  # 1. draw surrogate data
  S.chol <- S / sqrt(alpha)
  g <- rmvnorm(1, f, S.chol)
  
  # 2. compute implied latent variants
  C <- chol.fun(Sigma + S)
  R <- S - S %*% chol2inv(C) %*% S
  m <- R %*% S.inv %*% g
  
  L <- chol.fun(R)
  eta <- backsolve(L, f - m, transpose = TRUE)
  
  # 3. randomly center a bracket
  v <- runif(p, 0, sigma)
  theta.min <- theta - v
  theta.max <- theta.min + sigma
  
  # 4. draw u
  u <- runif(1)
  
  # 5. determine threshold
  log.y <- log.lik(u, f, C, g, theta, a, b, Y)
  
  # 6. draw proposal
  theta.s <- runif(p, theta.min, theta.max)
  
  # 7. compute function
  Sigma.s <- kern.fun(D, c(alpha, theta.s))
  C.s <- chol.fun(Sigma.s + S)
  R.s <- S - S %*% chol2inv(C.s) %*% S
  m.s <- R.s %*% S.inv %*% g
  
  L.s <- chol.fun(R.s)
  f.s <- crossprod(L.s, eta) + m.s
  
  log.yp <- log.lik(1, f.s, C.s, g, theta.s, a, b, Y)
  
  while (log.yp <= log.y) {
    for (pp in 1:p) {
      if (theta.s[pp] < theta[pp]) theta.min[pp] <- theta.s[pp] else theta.max[pp] <- theta.s[pp]
    }
    
    theta.s <- runif(p, theta.min, theta.max)
    
    Sigma.s <- kern.fun(D, c(alpha, theta.s))
    C.s <- chol.fun(Sigma.s + S)
    R.s <- S - S %*% chol2inv(C.s) %*% S
    m.s <- R.s %*% S.inv %*% g
    
    L.s <- chol.fun(R.s)
    f.s <- crossprod(L.s, eta) + m.s
    
    log.yp <- log.lik(1, f.s, C.s, g, theta.s, a, b, Y)
  }
  
  return(list("f" = f.s, "theta" = theta.s, "lp" = log.yp))
  
}

elliptical <- function(f0, C, Y) {
  # Samples f1 using elliptical slice sampler
  #
  # Args:
  #   f0       : current latent function
  #   C        : C'C = Sigma, f ~ GP(0, Sigma)
  #
  # Returns:
  #   New latent function f1
  
  log.L <- function(f) -sum(log1pe(-Y*f))
  
  n <- nrow(C)
  nu <- rmvnorm(1, rep(0, n), C)
  
  u <- runif(1)
  log.y <- log.L(f0) + log(u)
  
  theta <- runif(1, 0, 2*pi)
  theta.min <- theta - 2*pi; theta.max <- theta
  f1 <- f0 * cos(theta) + nu * sin(theta)
  
  log.L.f1 <- log.L(f1)
  
  while (log.L.f1 <= log.y) {
    if (theta < 0) theta.min <- theta else theta.max <- theta
    theta <- runif(1, theta.min, theta.max)
    f1 <- f0 * cos(theta) + nu * sin(theta)
    log.L.f1 <- log.L(f1)
  }
  
  return(list("f" = f1, "log.L" = log.L.f1))
  
}

gp.class <- function(dist, Y, split_index, a, b, ns = 1000, monitor = TRUE) {
  # Samples from p(f | X, Y)
  #
  # Args:
  #   dist        : distance array D s.t kernel K(D) with f ~ GP(0, K)
  #   Y           : Binary response vector
  #   split_index : Train/validate assignment
  #   a, b        : Hyperpriors for theta ~ InvGamma(a, b)
  #   ns          : number of samples
  #
  # Returns:
  #   Sample from latent posterior
  
  N.train <- sum(split_index == "train")
  train <- which(split_index == "train")
  D.train <- dist[train, train, , drop = FALSE]
  Y.train <- Y[train]
  
  # Initialize
  ## p(f | X, y)
  p <- length(a)
  theta <- matrix(nrow = p, ncol = ns)      # theta = (sigma, l_1, ..., l_p)
  theta[1, 1] <- 1                                      # sigma ~ IG(0, 0)
  theta[-1, 1] <- b[-1] / (a[-1] + 1)                   # ell ~ IG(a, b)
  f <- matrix(nrow = N.train, ncol = ns); f[, 1] <- 0   # prior mean
  lp.f <- numeric(ns); lp.f[1] <- 0
  
  for (t in 2:ns) {
    if (monitor && t %% (ns / 10) == 0) print(paste(t, "out of", ns))
    
    # jointly sample f and length scale(s)
    samp <- slice(theta = theta[-1, t - 1]
                  , f = f[, t - 1]
                  , alpha = theta[1, t - 1]
                  , Y = Y.train
                  , D = D.train
                  , a = a[-1], b = b[-1])
    
    f[, t] <- samp[["f"]]
    theta[-1, t] <- samp[["theta"]]
    lp.f[t] <- samp[["lp"]]
    
    # cheap updates of f
    K0 <- kern.fun(D.train, c(1, theta[-1, t])) # K0 = K / sigma
    K0.chol <- chol.fun(K0)
    K.chol <- K0.chol * sqrt(theta[1, t - 1])
    for (ii in 1:10) {
      slice <- elliptical(f[, t], K.chol, Y.train)
      f[, t] <- slice[["f"]]
    }
    lp.f[t] <- lp.f[t] + slice[["log.L"]]
    
    # sigma
    theta[1, t] <- rinvgamma(1, shape = a[1] + N.train / 2
                             , rate = crossprod(backsolve(K0.chol
                                                          , f[, t]
                                                          , transpose = TRUE)
                             ) / 2 + b[1])
    
    lp.f[t] <- lp.f[t] + dinvgamma(theta[1, t], shape = a[1] + N.train / 2
                                   , rate = crossprod(backsolve(K0.chol
                                                                , f[, t]
                                                                , transpose = TRUE)
                                   ) / 2 + b[1]
                                   , log = TRUE)
    K.chol <- K0.chol * sqrt(theta[1, t])
    
  }
  return(rbind(f, theta, t(as.matrix(lp.f))))
}

gp.class_pred <- function(gpc, D, split_index, p, burn = .2, thin = 10, avg = TRUE) {
  ns <- dim(gpc)[1]
  samp <- seq(burn*ns, ns, thin)
  
  N.train <- sum(split_index == "train")
  
  if (avg) {
    # posterior mean
    
    f.avg <- apply(gpc[samp, 1:N.train], 2, mean)
    
    # average hyperparameters
    theta.avg <- apply(gpc[samp, (N.train+1):(N.train+p+1)], 2, mean)
    K.avg <- kern.fun(D, replace(theta.avg, 1, 1)) # note that signal variance cancels
    fstar.avg <- crossprod(K.avg[split_index == "train"
                                 , split_index == "validate"]
                           , chol2inv(chol.fun(K.avg[split_index == "train"
                                                     , split_index == "train"]))) %*% f.avg
    y.avg <- sapply(fstar.avg, ilogit)
    
    return(y.avg)
    
  } else {
    # p(f* | X, y, x*)
    samp.fstar <- sapply(samp, function(ind) {
      theta.samp <- gpc[ind, (N.train+1):(N.train+p+1)]
      K <- kern.fun(D, replace(theta.samp, 1, 1))
      K.inv <- chol2inv(chol.fun(K[split_index == "train", split_index == "train"]))
      fstar <- crossprod(K[split_index == "train", split_index == "validate"]
                         , K.inv) %*% gpc[ind, 1:N.train]
    })
    # p(y* = 1 | f*, y, x*)
    if (is.matrix(samp.fstar)) {
      ystar <- apply(samp.fstar, 2, function(i) sapply(i, ilogit))
      y.map <- rowMeans(ystar)
    } else {
      ystar <- sapply(samp.fstar, ilogit)
      y.map <- mean(ystar)
    }
    
    return(y.map)
  }
}

mcmc_array <- function (ns, nchains, params) {
  nparams <- length(params)
  array(dim = c(ns, nchains, nparams),
        dimnames = list(iterations = NULL,
                        chains = paste0("chain", 1:nchains),
                        parameters = params))
}

mcc <- function(y, y.pred) {
  
  y.pred <- factor(y.pred, levels = c(-1, 1)) # in case only one class predicted
  
  cm <- table(y, y.pred)
  
  TP <- cm["1", "1"]
  TN <- cm["-1", "-1"]
  FP <- cm["-1", "1"]
  FN <- cm["1", "-1"]
  
  num <- TP*TN - FP*FN
  den <- (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)
  den <- ifelse(den == 0, 1, den)
  
  return(num/sqrt(den)) 
}

roc <- function(labels, scores, plot = FALSE) {
  
  if (min(labels) == -1) labels <- (labels + 1) / 2
  
  Labels <- labels[order(scores, decreasing = TRUE)]
  FPR <- cumsum(!Labels)/sum(!Labels)
  TPR <- cumsum(Labels)/sum(Labels)
  
  if (plot) {
    plot(TPR ~ FPR, type = "o")
    abline(c(0, 1))
  }
  
  N <- length(labels)
  AUC <- sum((FPR[2:N] - FPR[1:N-1]) * TPR[2:N])
  
  return(AUC)
}

f1 <- function(y, y.pred) {
  
  y.pred <- factor(y.pred, levels = c(-1, 1)) # in case only one class predicted
  
  cm <- table(y, y.pred)
  
  TP <- cm["1", "1"]
  TN <- cm["-1", "-1"]
  FP <- cm["-1", "1"]
  FN <- cm["1", "-1"]
  
  precision <- TP / (TP + FP)
  recall <- TP / (TP + FN)
  
  F1 <- 2 * (precision * recall) / (precision + recall)
  
  return(F1) 
}