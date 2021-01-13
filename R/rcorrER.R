rcorrER <- function(m = 1, G, s) {
  # random samples from correlated ER graph
  #
  # m : number of samples
  # G : true underlying graph
  # s : correlation
  
  library(Matrix)
  
  if (is.igraph(G)) G <- get.adjacency(G)
  
  n <- nrow(G)
  true.edges <- which(G == 1)
  
  G_i <- replicate(m, expr = { # draw m samples
    
    # Add true edges from G but thin wp s
    g <- G
    g[true.edges] <- sample(0:1, length(true.edges)
                            , replace = TRUE
                            , prob = c(s, 1 - s))
    
    # Make simple and undirected
    diag(g) <- 0
    g <- forceSymmetric(g, uplo = "L")
    
    return(g)
  })
  
  if (m == 1) G_i <- G_i[[1]]
  
  return(G_i)
}
