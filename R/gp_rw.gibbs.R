gp_rw.class <- function(K, Y, split_index, ns = 1000, monitor = TRUE) {
  # Samples from p(f | X, Y)
  #
  # Args:
  #   K           : RW kernel K with f ~ GP(0, K)
  #   Y           : Binary response vector
  #   split_index : Train/validate assignment
  #   ns          : number of samples
  #
  # Returns:
  #   Sample from latent posterior
  
  N.train <- sum(split_index == "train")
  train <- which(split_index == "train")
  Y.train <- Y[train]
  
  # Initialize
  ## p(f | X, y)
  K.chol <- chol.fun(K[split_index == "train", split_index == "train"])
  f <- matrix(nrow = N.train, ncol = ns); f[, 1] <- 0   # prior mean
  lp.f <- numeric(ns); lp.f[1] <- 0
  
  for (t in 2:ns) {
    if (monitor && t %% (ns / 10) == 0) print(paste(t, "out of", ns))
    
    slice <- elliptical(f[, t-1], K.chol, Y.train)
    f[, t] <- slice[["f"]]
    lp.f[t] <- lp.f[t] + slice[["log.L"]]
    
    lp.f[t] <- lp.f[t]
  }
  return(rbind(f, t(as.matrix(lp.f))))
}

gp_rw.class_pred <- function(gpc, K, split_index, burn = .2, thin = 10, avg = TRUE) {
  ns <- dim(gpc)[1]
  samp <- seq(burn*ns, ns, thin)
  
  N.train <- sum(split_index == "train")
  
  if (avg) {
    # posterior mean
    f.avg <- apply(gpc[samp, 1:N.train], 2, mean)
    fstar.avg <- crossprod(K[split_index == "train"
                             , split_index == "validate"]
                           , chol2inv(chol.fun(K[split_index == "train"
                                                 , split_index == "train"]))) %*% f.avg
    y.avg <- sapply(fstar.avg, ilogit)
    
    return(y.avg)
    
  } else {
    # p(f* | X, y, x*)
    K.inv <- chol2inv(chol.fun(K[split_index == "train", split_index == "train"]))
    samp.fstar <- sapply(samp, function(ind) {
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