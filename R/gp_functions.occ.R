gp.occ <- function(gpc, burn = .2, thin = 10) {
  # Computes OCC scores from GP classifer
  #
  # Args:
  #   gpc         : gpc.class object
  #   burn        : how long is burn-in
  #   thin        : how often to thin
  #
  # Returns:
  #   OCC scores

  ns <- dim(gpc)[1]
  samp <- seq(burn*ns, ns, thin)
  
  # p(f* | X, y, x*)
  samp.fstar <- sapply(samp, function(ind) {
    theta.samp <- gpc[ind, (N.train+1):(N.train+p+1)]
    K <- kern.fun(D, replace(theta.samp, 1, 1))
    K.inv <- chol2inv(chol.fun(K[split_index == "train", split_index == "train"]))
    fstar <- crossprod(K[split_index == "train", split_index == "validate"]
                       , K.inv) %*% gpc[ind, 1:N.train]
  })
  
  # p(y* = 1 | f*, y, x*)
  ystar <- apply(samp.fstar, 2, function(i) sapply(i, ilogit))
  y.map <- rowMeans(ystar)
  
  # average hyperparameters
  f.avg <- apply(gpc[samp, 1:N.train], 2, mean)
  
  theta.avg <- apply(gpc[samp, (N.train+1):(N.train+p+1)], 2, mean)
  K.avg <- kern.fun(D, replace(theta.avg, 1, 1)) # note that signal variance cancels
  fstar.avg <- crossprod(K.avg[split_index == "train"
                               , split_index == "validate"]
                         , chol2inv(chol.fun(K.avg[split_index == "train"
                                                   , split_index == "train"]))) %*% f.avg
  y.avg <- sapply(fstar.avg, ilogit)
  
  # variance
  sigma_star <- apply(ystar, 1, sd)
  
  return(list(mu = fstar.avg, pi = y.avg, sigma = sigma_star, H = fstar.avg / sqrt(sigma_star)))
  
}

plot.occ <- function(occ.scores) {
  require(ggplot2)
  require(gridExtra)
  
  p1 <- ggplot(data.frame(mu = sort(occ.scores$mu)
                          , ord = 1:length(occ.scores$mu)
                          , Class = as.factor(Y[split_index == "validate"][order(occ.scores$mu)]))
               , aes(x = ord, y = mu, color = Class)) +
    geom_point(aes(shape = Class)) +
    geom_hline(aes(yintercept = mean(sort(occ.scores$mu)[c(which.max(diff(sort(occ.scores$mu)))
                                                           , which.max(diff(sort(occ.scores$mu))) + 1)]))
               , linetype = "dashed") +
    scale_shape_manual(values = c(1, 3)) + 
    scale_color_manual(values = c("black", "red")) +
    theme_bw() +
    ylab(expression(tilde(mu))) +
    theme(axis.title.x = element_blank()
          , axis.text.x = element_blank()
          , axis.ticks.x = element_blank()
          , axis.title.y = element_text(size = 20)
          , axis.ticks.y = element_blank()
          , axis.text.y = element_blank()
          , legend.position = "none")
  
  p2 <- ggplot(data.frame(pi = sort(occ.scores$pi)
                          , ord = 1:length(occ.scores$pi)
                          , Class = as.factor(Y[split_index == "validate"][order(occ.scores$pi)]))
               , aes(x = ord, y = pi, color = Class)) +
    geom_point(aes(shape = Class)) +
    geom_hline(aes(yintercept = mean(sort(occ.scores$pi)[c(which.max(diff(sort(occ.scores$pi)))
                                                           , which.max(diff(sort(occ.scores$pi))) + 1)]))
               , linetype = "dashed") +
    scale_shape_manual(values = c(1, 3)) + 
    scale_color_manual(values = c("black", "red")) +
    theme_bw() +
    ylab(expression(tilde(pi))) +
    theme(axis.title.x = element_blank()
          , axis.text.x = element_blank()
          , axis.ticks.x = element_blank()
          , axis.title.y = element_text(size = 20)
          , axis.ticks.y = element_blank()
          , axis.text.y = element_blank()
          , legend.position = "none")
  
  p3 <- ggplot(data.frame(sigma = sort(occ.scores$sigma)
                          , ord = 1:length(occ.scores$sigma)
                          , Class = as.factor(Y[split_index == "validate"][order(occ.scores$sigma)]))
               , aes(x = ord, y = sigma, color = Class)) +
    geom_point(aes(shape = Class)) +
    geom_hline(aes(yintercept = mean(sort(occ.scores$sigma)[c(which.max(diff(sort(occ.scores$sigma)))
                                                              , which.max(diff(sort(occ.scores$sigma))) + 1)]))
               , linetype = "dashed") +
    scale_shape_manual(values = c(1, 3)) + 
    scale_color_manual(values = c("black", "red")) +
    theme_bw() +
    ylab(expression(tilde(sigma))) +
    theme(axis.title.x = element_blank()
          , axis.text.x = element_blank()
          , axis.ticks.x = element_blank()
          , axis.title.y = element_text(size = 20)
          , axis.ticks.y = element_blank()
          , axis.text.y = element_blank()
          , legend.position = "none")
  
  p4 <- ggplot(data.frame(H = sort(occ.scores$H)
                          , ord = 1:length(occ.scores$H)
                          , Class = as.factor(Y[split_index == "validate"][order(occ.scores$H)]))
               , aes(x = ord, y = H, color = Class)) +
    geom_point(aes(shape = Class)) +
    geom_hline(aes(yintercept = mean(sort(occ.scores$H)[c(which.max(diff(sort(occ.scores$H)))
                                                          , which.max(diff(sort(occ.scores$H))) + 1)]))
               , linetype = "dashed") +
    scale_shape_manual(values = c(1, 3)) + 
    scale_color_manual(values = c("black", "red")) +
    theme_bw() +
    ylab(expression(H)) +
    theme(axis.title.x = element_blank()
          , axis.text.x = element_blank()
          , axis.ticks.x = element_blank()
          , axis.title.y = element_text(size = 20)
          , axis.ticks.y = element_blank()
          , axis.text.y = element_blank()
          , legend.position = "none")
  
  grid.arrange(p1, p2, p3, p4, ncol = 2)
  
}