# Libraries ----
library(igraph)
library(invgamma) # use this implementation of rinvgamma
# library(actuar) # implementation of rinvgamma switches rate and scale
library(abind)
library(Matrix)

# Functions ----
ham_dis <- function(G1, G2) {
  # computes modified Hamming distance
  
  n <- nrow(G1)
  sum((G1 - G2)^2) / (n * (n-1))       # square (equivalent for binary)
}

kern.dis <- function(G, X = NULL) {
  # computes distance matrix (array) given G (and X)
  # NB: distance matrix for X is squared-Euclidean distance
  
  require(parallel)
  
  p.G <- ifelse(missing(G), 0, length(G))
  p.X <- ifelse(is.null(ncol(X)), 0, ncol(X))
  
  if (!missing(G)) {
    N <- length(G[[1]])
    D <- array(0, dim = c(N, N, p.G + p.X))
    
    for (pp in 1:p.G) {
      D.ham <- matrix(0, N, N)
      d <- unlist(mclapply(mc.cores = detectCores() - 1
                           , 1:(N - 1)
                           , function(i) sapply((i+1):N
                                                , function(j)
                                                  ham_dis(G[[pp]][[i]]
                                                          , G[[pp]][[j]])
                           )))
      
      D.ham[lower.tri(D.ham, diag = F)] <- d
      D.ham[upper.tri(D.ham, diag = F)] <- t(D.ham)[upper.tri(D.ham, diag = F)]
      D[, , pp] <- D.ham
    }
  }
  
  if (!is.null(X)) {
    N <- nrow(X); p <- ncol(X)
    
    if (missing(G)) {
      D <- array(0, dim = c(N, N, p.X))
    }
    
    for (pp in 1:p.X) {
      D.pp <- matrix(0, N, N)
      d <- unlist(mclapply(mc.cores = detectCores() - 1
                           , 1:(N - 1)
                           , function(i) sapply((i+1):N
                                                , function(j)
                                                  (X[i, pp] - X[j, pp])^2
                           )))
      
      D.pp[lower.tri(D.pp, diag = F)] <- d
      D.pp[upper.tri(D.pp, diag = F)] <- t(D.pp)[upper.tri(D.pp, diag = F)]
      D[, , p.G + pp] <- D.pp
    }
  }
  
  return(D)
}

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

chol.fun <- function(x, jit = 1e-2) {
  # compute Cholesky function with added jitter
  
  n <- nrow(x)
  L <- chol(x + diag(n) * jit^2)
  
  return(L)
}

rmvnorm <- function(n = 1, mu, C){
  # samples from MVN given Cholesky of covariance

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

# Survival Functions ----

kern.dis_star <- function(old, new, square = FALSE) {
  # computes distance matrix between old and new (Y values)
  # NB: default is not to return full square matrix
  
  n <- length(old)
  m <- length(new)
  D <- as.matrix(dist(c(old, new)))
  D <- D^2
  
  if (square) {
    return(D)
  } else {
    return(as.matrix(D[1:n, (n+1):(n+m)]))
  }
  
}

kern.dis_expand_grid <- function(D, G, grid, ind) {
  # expands distance matrix (array) for tracking covariates
  
  n <- dim(D)[1]
  p <- dim(D)[3]
  n_i <- length(grid)
  n_G_i <- unlist(lapply(G, length))
  n_G <- sum(n_G_i)
  n_G_cumsum <- c(0, cumsum(n_G_i))
  
  D.rep_vec <- array(0, c(1, n+n_G, p))
  D.rep_vec[1, 1:n, ] <- D[ind, , , drop = FALSE]
  
  if (n_G > 0) {
    for (i in 1:n) {
      if (n_G_i[i] == 0) {
        next
      }
      for (pp in 1:p) {
        D.rep_vec[1, (n+n_G_cumsum[i]+1):(n+n_G_cumsum[i+1]), pp] <- D[ind, i, pp]
      }
    }
  }
  
  D.rep <- array(0, c(n+n_G, n_i, p))
  for (pp in 1:p) {
    D.rep[, , pp] <- t(matrix(D.rep_vec[, , pp], nrow = n_i, ncol = n+n_G, byrow = TRUE))
  }
  
  return(D.rep)
}

kern.dis_expand_A <- function(D, G, A, ind) {
  # expands distance matrix (array) for tracking covariates
  
  n <- dim(D)[1]
  p <- dim(D)[3]
  n_i <- length(A)
  n_G_i <- unlist(lapply(G, length))
  n_G <- sum(n_G_i)
  n_G_cumsum <- c(0, cumsum(n_G_i))
  
  D.rep_vec <- array(0, c(1, n+n_G, p))
  D.rep_vec[1, 1:n, ] <- D[ind, , , drop = FALSE]
  
  if (n_G > 0) {
    for (i in 1:n) {
      if (n_G_i[i] == 0) {
        next
      }
      for (pp in 1:p) {
        D.rep_vec[1, (n+n_G_cumsum[i]+1):(n+n_G_cumsum[i+1]), pp] <- D[ind, i, pp]
      }
    }
  }
  
  D.rep <- array(0, c(n+n_G, n_i, p))
  for (pp in 1:p) {
    D.rep[, , pp] <- t(matrix(D.rep_vec[, , pp], nrow = n_i, ncol = n+n_G, byrow = TRUE))
  }
  
  return(D.rep)
}

kern.dis_expand <- function(D, G_i) {
  # expands distance matrix (array) for tracking covariates
  
  n <- dim(D)[1]
  p <- dim(D)[3]
  
  n.G_i <- sapply(G_i, length)
  n.G_i.sum <- c(0, cumsum(n.G_i))
  n.G <- sum(n.G_i)
  
  D.rep <- array(0, c(n+n.G, n+n.G, p))
  D.rep[1:n, 1:n,] <- D
  for (i in 1:(n-1)) {
    if (n.G_i[i] == 0) {
      next
    }
    for (j in (i+1):n) {
      if (n.G_i[j] == 0) {
        next
      }
      for (pp in 1:p) {
        D.rep[(n+n.G_i.sum[i]+1):(n+n.G_i.sum[i+1])
              , (n+n.G_i.sum[j]+1):(n+n.G_i.sum[j+1])
              , pp] <- matrix(D[i, j, pp], n.G_i[i], n.G_i[j])
      }
    }
  }
  for (j in 1:n) {
    if (n.G_i[j] == 0) {
      next
    }
    for (pp in 1:p) {
      D.rep[1:n
            , (n+n.G_i.sum[j]+1):(n+n.G_i.sum[j+1])
            , pp] <- matrix(D[, j, pp], n, n.G_i[j])
    }
  }
  for (pp in 1:p) {
    D.rep[,,pp] <- as.matrix(forceSymmetric(D.rep[,,pp]))
  }
  
  return(D.rep)
}

trap_rule <- function(grid, grid.values) {
  # integration using trapezoidal rule

  n <- length(grid.values)
  delta <- grid[2] - grid[1]
  sum.pairs <- grid.values[-n] + grid.values[-1]
  return(delta/2 * cumsum(sum.pairs))
}

comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

slice <- function(theta, f, alpha, n.G, D, a, b, sigma = 10) {
  # Implementation from:
  # "Slice sampling covariance hyperparameters of latent Gaussian models"
  #
  # Jointly samples hyperparameters theta and latent function f
  #
  # Args:
  #   theta       : kernel hyperparameters
  #   f           : current latent position
  #   alpha       : auxiliary noise
  #   n.G         : |G|, i.e number of rejected points
  #   D           : distance array D s.t kernel K(D) with f ~ GP(0, K)
  #   a, b        : Hyperpriors for theta ~ Inv-Gamma(a, b)
  #   sigma       : width of slice
  #
  # Returns:
  #   New theta and f
  
  log.lik <- function(u, f, C, g, theta, a, b) {
    log(u) +                                                 # u
      -sum(log1pe(-f)) - sum(f[(n-n.G+1):n]) +               # L(f)
      -.5 * (crossprod(backsolve(C, g, transpose = TRUE))) + # N(g; 0, Sigma + S)
      sum(-a * theta - b / exp(theta))                       # p(theta)
  }
  
  n <- length(f)
  p <- length(theta)
  
  Sigma <- kern.fun(D, c(alpha, theta))
  S <- diag(n) * alpha
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
  log.y <- log.lik(u, f, C, g, theta, a, b)
  
  # 6. draw proposal
  theta.s <- runif(p, theta.min, theta.max)
  
  # 7. compute function
  Sigma.s <- kern.fun(D, c(alpha, theta.s))
  C.s <- chol.fun(Sigma.s + S)
  R.s <- S - S %*% chol2inv(C.s) %*% S
  m.s <- R.s %*% S.inv %*% g
  
  L.s <- chol.fun(R.s)
  f.s <- crossprod(L.s, eta) + m.s
  
  log.yp <- log.lik(1, f.s, C.s, g, theta.s, a, b)
  
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
    
    log.yp <- log.lik(1, f.s, C.s, g, theta.s, a, b)
  }
  
  return(list("f" = f.s, "theta" = theta.s, "lp" = log.yp))
  
}

elliptical <- function(f0, n.G, C) {
  # Implementation from:
  # "Elliptical slice sampling"
  #
  # Samples f1 using elliptical slice sampler
  #
  # Args:
  #   f0       : current latent function
  #   n.G      : |G|
  #   C        : C'C = Sigma, f ~ GP(0, Sigma)
  #
  # Returns:
  #   New latent function f1
  
  n <- nrow(C)
  
  log.L <- function(f) -sum(log1pe(-f)) - sum(f[(n-n.G+1):n])
  
  nu <- rmvnorm(1, rep(0, n), C)
  
  u <- runif(1)
  log.y <- log.L(f0) + log(u)
  
  theta <- runif(1, 0, 2*pi)
  theta.min <- theta - 2*pi; theta.max <- theta
  f1 <- f0 * cos(theta) + nu * sin(theta)
  
  while (log.L(f1) <= log.y) {
    if (theta < 0) theta.min <- theta else theta.max <- theta
    theta <- runif(1, theta.min, theta.max)
    f1 <- f0 * cos(theta) + nu * sin(theta)
  }
  
  return(f1)
  
}

gp.survival <- function(dist, Y, a, b, N = 100
                        , burn = .2, thin = 10, grid.length = 100) {
  # Algorithm 2 from paper
  #
  # Samples from GP Survival posterior 
  # with exponential baseline hazard function
  #
  # Args:
  #   dist        : distance array D s.t kernel K(D) with f ~ GP(0, K)
  #   Y           : survival times
  #   a, b        : Hyperpriors for Omega and hyperparameters
  #   N           : number of samples
  #   burn        : how long is burn-in
  #   thin        : how often to thin
  #   grid.length : grid length for trapezoidal rule
  #
  # Returns:
  #   Chains for hyperparameters and survival
  
  # Initialize
  
  require(foreach)
  require(doParallel)
  
  D.mat <- dist
  n <- dim(D.mat)[1]
  p <- dim(D.mat)[3]
  
  hyperparams <- matrix(nrow = p+1, ncol = N)
  if (p == 2) {
    rownames(hyperparams) <- c("sigma", "ell.net", "ell.time")
  } else { # non-network covariates
    rownames(hyperparams) <- c("sigma", "ell.net"
                               , paste0("ell.X", 1:(p-2))
                               , "ell.time") # must be last
  }
  
  hyperparams["sigma", 1] <- 1                      # sigma ~ flat
  hyperparams[-1, 1] <- b[-(1:2)] / (a[-(1:2)] + 1) # ell ~ IG(a, b)
  
  Omega <- vector(length = N)
  Omega[1] <- 1 / mean(Y)
  Y.sum <- sum(Y)
  
  n_i <- vector(length = n)
  A_i <- vector("list", n)
  
  l <- vector("list", N)
  K <- kern.fun(D.mat, hyperparams[, 1])
  K.chol <- chol.fun(K)
  l[[1]] <- rep(0, n)
  l.A <- vector("list", n)
  U_i <- vector("list", n)
  G <- NULL
  
  delta <- max(Y) * 1.025 / (grid.length - 1)
  grid <- seq(0, max(Y) * 1.025, delta)
  store <- round(seq(burn*N, N, thin))
  surv <- array(dim = c(grid.length, n, length(store)))
  K.chol <- chol.fun(K)
  K.inv <- chol2inv(K.chol)
  D.grid <- kern.dis_star(Y, grid)
  
  cores <- detectCores()
  cl <- makeCluster(cores[1] - 1)
  registerDoParallel(cl)
  
  for (q in 2:N) {
    if (q %% floor(N / 10) == 0) print(paste(q, "out of", N))
    
    # lines 2-6 : mapping homogeneous Poisson process into a 
    #             non-homogeneous with the appropriate intensity
    
    lambda_0 <- Omega[q - 1]                # constant baseline hazard : lambda_0 ~ Exp(Omega)
    
    samp <- foreach(i = 1:n, .combine = "comb", .multicombine = TRUE
                    , .init = list(list(), list()), .packages = "abind"
                    , .export = c("kern.dis_expand_A", "kern.dis_star", "kern.fun"
                                  , "rmvnorm", "chol.fun", "log1pe")
    ) %dopar% {
      
      Lambda_0 <- lambda_0*Y[i]             # cumulative hazard function
      n_i <- rpois(1, Lambda_0)
      A.tilde <- runif(n_i, max = Lambda_0)
      A_i <- A.tilde / lambda_0        # Lambda^{-1} (A.tilde)
      
      # line 7 : sample l(A) | l(G U T), lambda_0
      if (n_i == 0) { # no A_i generated from Poisson process
        return(list(c(), c()))
        
      } else {
        
        D.A_rep <- kern.dis_expand_A(dist[, , -p, drop = FALSE], G, A_i, i)
        D.A <- kern.dis_star(c(Y, unlist(G)), A_i)
        D.A_mat <- abind(D.A_rep, D.A)
        Ks <- kern.fun(D.A_mat, hyperparams[, q-1])
        
        K.A <- kern.fun(array(as.matrix(dist(A_i))^2, dim = c(n_i, n_i, 1))
                        , hyperparams[c("sigma", "ell.time"), q-1]) +
          (p-1) * hyperparams["sigma", q-1] # all within same covariates
        K.A <- K.A - crossprod(backsolve(K.chol, Ks, transpose = TRUE))
        
        l.A <- rmvnorm(mu = crossprod(Ks, chol2inv(K.chol)) %*% l[[q-1]]
                       , C = chol.fun(K.A))
        
        # lines 8-11 : update G
        U_i <- runif(n_i)
        # G_i <- A_i[U_i < 1 - sapply(l.A, ilogit)]
        G_i <- A_i[log(U_i) > -log1pe(-l.A)] # u ~ 1 - u
        l.G <- l.A[log(U_i) > -log1pe(-l.A)]
        
        return(list(G_i, l.G))
      }
      
    } # end for loop
    
    G_i <- samp[[1]]
    l.G <- unlist(samp[[2]])
    
    # line 12 : update parameters of hazard function
    # Omega
    n.G <- sum(sapply(G_i, length))
    Omega[q] <- rgamma(1
                       , shape = a[1] + n + n.G
                       , rate = b[1] + Y.sum)
    
    # line 13 : update l(G U T) and hyperparameters
    # l and ell jointly
    G <- G_i
    D.net_rep <- kern.dis_expand(dist[, , -p, drop = FALSE], G)
    D.G <- kern.dis_star(Y, unlist(G), square = TRUE)
    D.mat <- abind(D.net_rep, D.G)
    
    samp <- slice(theta = hyperparams[-1, q - 1]
                  , f = c(l[[q-1]][1:n], l.G)
                  , alpha = hyperparams["sigma", q - 1]
                  , n.G = n.G
                  , D = D.mat
                  , a = a[-(1:2)], b = b[-(1:2)])
    
    l[[q]] <- samp[["f"]]
    hyperparams[-1, q] <- samp[["theta"]]
    
    # cheap updates of l
    K0 <- kern.fun(D.mat, c(1, hyperparams[-1, q])) # K0 = K / sigma
    K0.chol <- chol.fun(K0)
    K.chol <- K0.chol * sqrt(hyperparams["sigma", q - 1])
    for (ii in 1:10) {
      l[[q]] <- elliptical(l[[q]], n.G, K.chol)
    }
    
    # sigma
    hyperparams["sigma", q] <- rinvgamma(1, shape = a[2] + (n + n.G) / 2
                                         , rate = crossprod(backsolve(K0.chol
                                                                      , l[[q]]
                                                                      , transpose = TRUE)
                                         ) / 2 + b[2])
    
    K.chol <- K0.chol * sqrt(hyperparams["sigma", q])
    
    if (q %in% store) {
      # l(grid)
      K.inv <- chol2inv(K.chol)
      D.grid <- kern.dis_star(c(Y, unlist(G)), grid)
      for (j in 1:n) {
        D.grid_rep <- kern.dis_expand_grid(D[, , -p, drop = FALSE], G, grid, j)
        D.grid_mat <- abind(D.grid_rep, D.grid)
        Ks <- kern.fun(D.grid_mat, hyperparams[, q])
        surv[, j, which(store == q)] <- crossprod(Ks, K.inv) %*% l[[q]]
      } 
    }
    
  }
  
  on.exit(stopCluster(cl))
  
  return(list(l = l
              , Omega = Omega
              , hyperparams = hyperparams
              , surv = surv))
}

plot.surv <- function(gp.surv, Y, groups, burn = .2, thin = 10, km = FALSE) {
  
  require(ggplot2)
  
  N <- length(Y)
  ns <- length(gp.surv$l)
  grid.length <- dim(gp.surv$surv)[1]
  grid <- seq(0, max(Y)*1.025, max(Y)*1.025 / (grid.length - 1))
  store <- seq(burn*ns, ns, thin)
  
  hazard.post <- sapply(1:N, function(i)
    rowMeans(matrix(gp.surv$Omega[store], nrow = grid.length, ncol = length(store), byrow = TRUE)
             * sapply(gp.surv$surv[, i, ], ilogit)))
  
  # Individual Survival Surfaces
  int.trap <- trap_rule(grid, hazard.post[, 1])
  surv.j <- exp(-c(0, int.trap))
  df <- data.frame(surv = surv.j, grid, Class = as.factor(groups[1]), subject = as.factor(1))
  for (j in 2:N) {
    int.trap <- trap_rule(grid, hazard.post[, j])
    surv.j <- exp(-c(0, int.trap))
    df <- rbind(df, data.frame(surv = surv.j, grid, Class = as.factor(groups[j]), subject = as.factor(j)))
  }
  
  p1 <- ggplot(df, aes(x = grid, y = surv, col = Class)) +
    geom_line(aes(group = subject), alpha = .1) +
    theme_bw() +
    geom_rug(data = data.frame(Y, groups), aes(x = Y, col = as.factor(groups))
             , sides = "b", inherit.aes = FALSE) +
    labs(x = "t", y = "S(t)") +
    scale_colour_manual(values = c("#F8766D", "#00BFC4")
                        , labels = expression(f[0], f[1]))
  
  # Grouped Survival Curves
  hazard.post.0 <- rowMeans(hazard.post[, groups == 0])
  int.trap <- trap_rule(grid, hazard.post.0)
  surv.0 <- exp(-c(0, int.trap))
  hazard.post.1 <- rowMeans(hazard.post[, groups == 1])
  int.trap <- trap_rule(grid, hazard.post.1)
  surv.1 <- exp(-c(0, int.trap))
  df <- data.frame(surv = c(surv.0, surv.1), grid
                   , Class = as.factor(c(rep(0, length(grid))
                                         , rep(1, length(grid)))))
  
  p1 <- p1 + geom_line(data = df, aes(x = grid, y = surv, col = Class), size = 1.5)
  
  if (km) {
    require(survival)
    
    km <- survfit(Surv(Y) ~ groups)
    df <- data.frame(surv = km$surv
                     , grid = km$time
                     , group = as.factor(c(rep(0, which.max(diff(km$surv)) - 1), rep(1, length(km$surv) - which.max(diff(km$surv)) + 1))))
    p1 <- p1 + geom_step(data = df, aes(x = grid, y = surv, group = group, color = group)
                         , inherit.aes = FALSE, size = 1, linetype = "longdash")
  }
  
  p1
}