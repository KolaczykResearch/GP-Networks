dist.eigen <- function(G, X = NULL, normalized = FALSE, signed = FALSE) {
  
  p.G <- ifelse(missing(G), 0, length(G))
  p.X <- ifelse(is.null(ncol(X)), 0, ncol(X))
  
  if (!missing(G)) {
    N <- length(G[[1]])
    D <- array(0, dim = c(N, N, p.G + p.X))
    n <- nrow(G[[1]][[1]])
    
    for (pp in 1:p.G) {
      
      if (!signed) {
        if (normalized) { # L = I - D^-1/2 A D^-1/2
          
          n <- nrow(G[[1]][[1]])
          L <- lapply(G[[pp]], function(A) {
            
            
            D <- diag(1 / sqrt(apply(A, 1, sum)))
            D[is.infinite(D)] <- 0
            I = diag(n)
            L <- I - D %*% A %*% D
          })
          
        } else { # L = D - A
          
          L <- lapply(G[[pp]], function(A) diag(apply(A, 1, sum)) - A)
        }
      } else { # signed
        
        if (normalized) { # L = I - D^-1/2 A D^-1/2
          
          n <- nrow(G[[1]][[1]])
          L <- lapply(G[[pp]], function(A) {
            
            
            D <- diag(1 / sqrt(apply(A, 1, function(x) sum(abs(x)))))
            D[is.infinite(D)] <- 0
            I = diag(n)
            L <- I - D %*% A %*% D
          })
          
        } else { # L = |D| - A
          
          L <- lapply(G[[pp]], function(A) diag(apply(A, 1, function(x) sum(abs(x)))) - A)  
        }
        
      }
      
      lambda <- sapply(L, function(x) sort(eigen(x, only.values = TRUE
                                                 , symmetric = TRUE)$values))
      D.eigen <- as.matrix(dist(t(lambda)))
      D[, , pp] <- D.eigen
    }
  }
  
  if (!is.null(X)) {
    N <- nrow(X); p <- ncol(X)
    
    if (missing(G)) {
      D <- array(0, dim = c(N, N, p.X))
    }
    
    for (pp in 1:p.X) {
      D[, , p.G + pp] <- as.matrix(dist(X[, pp])^2)
    }
  }
  
  return(D)
}

dist.frobenius <- function(G, X = NULL) {
  
  p.G <- ifelse(missing(G), 0, length(G))
  p.X <- ifelse(is.null(ncol(X)), 0, ncol(X))
  
  if (!missing(G)) {
    N <- length(G[[1]])
    D <- array(0, dim = c(N, N, p.G + p.X))
    n <- nrow(G[[1]][[1]])
    
    for (pp in 1:p.G) {
      V <- sapply(G[[pp]], as.vector)
      D.frobenius <- as.matrix(dist(t(V))^2) / (n*(n-1))
      D[, , pp] <- D.frobenius
    }
  }
  
  if (!is.null(X)) {
    N <- nrow(X); p <- ncol(X)
    
    if (missing(G)) {
      D <- array(0, dim = c(N, N, p.X))
    }
    
    for (pp in 1:p.X) {
      D[, , p.G + pp] <- as.matrix(dist(X[, pp])^2)
    }
  }
  
  return(D)
}