# Libraries ----
library(SpiecEasi)
library(ggplot2)
library(grid)
library(dplyr)
library(foreach)
library(doParallel)

# Data ----
# from : http://statweb.stanford.edu/~susan/papers/PNASRR.html
load("/Data/PregnancyClosed15.Rdata")

df <- data.frame(SubjectID = PS@sam_data@.Data[[5]]
                 , GDDel = PS@sam_data@.Data[[19]]
                 , GWDel = PS@sam_data@.Data[[20]]
                 , Outcome = PS@sam_data@.Data[[62]]
                 , History = PS@sam_data@.Data[[12]])
df <- unique(df)
Outcomes <- df[, c("SubjectID", "Outcome", "GWDel", "GDDel")]
History <- df[, c("SubjectID", "History")]

# Figure 1 ----
df.hist <- merge(Outcomes, History)
df.hist$Outcome <- factor(df.hist$Outcome
                          , levels = c("VeryPreterm", "Preterm", "Marginal", "Term"))
df.hist$History <- factor(df.hist$History, levels = c("1", "0"))
  
p <- ggplot(df.hist, aes(x = GWDel, fill = History)) +
  geom_histogram(binwidth = 1, color = "black") +
  labs(x = "Gestational Weeks", y = "Count") +
  scale_fill_manual(values = c("#F8766D", "#00BFC4")
                    , name = ""
                    , labels = c("History", "No History")) +
  scale_x_continuous(breaks = seq(22, 42, 2)) +
  geom_segment(aes(x = 37.5, y = 0, xend = 37.5, yend = Inf), linetype = "longdash") +
  theme_bw() +
  theme(plot.margin = unit(c(1, 1, 2, 1), "lines"))

Text1 = textGrob("Preterm")
Text2 = textGrob("Term")

p1 = p + 
  annotation_custom(grob = linesGrob(), xmin = 21, xmax = 21, ymin = -4.5, ymax = -3.75) +
  annotation_custom(grob = linesGrob(), xmin = 21, xmax = 37.5, ymin = -4.5, ymax = -4.5) +
  annotation_custom(grob = linesGrob(), xmin = 37.5, xmax = 37.5, ymin = -4.5, ymax = -3.75) +
  annotation_custom(grob = Text1,  xmin = 21, xmax = 37.5, ymin = -5.5, ymax = -5.5)

p1 = p1 + 
  annotation_custom(grob = linesGrob(), xmin = 37.5, xmax = 43, ymin = -4.5, ymax = -4.5) +
  annotation_custom(grob = linesGrob(), xmin = 43, xmax = 43, ymin = -4.5, ymax = -3.75) +
  annotation_custom(grob = Text2,   xmin = 37.5, xmax = 43, ymin = -5.5, ymax = -5.5) 

gt <- ggplot_gtable(ggplot_build(p1))
gt$layout$clip[gt$layout$name=="panel"] <- "off"
grid.draw(gt)

# Network Construction ----
weeks <- 1:30
df <- subset_samples(PS, GWColl %in% weeks)

otu <- otu_table(df)@.Data
otu.present <- otu; otu.present[otu > 0] <- 1
taxa_rare <- which(colSums(otu.present) < 500 | colSums(otu) < 1500)

samples <- data.frame(SubjectID = df@sam_data@.Data[[5]]
                      , GDColl = df@sam_data@.Data[[8]]
                      , GWColl = df@sam_data@.Data[[9]]
                      , Preg = df@sam_data@.Data[[56]]
                      , BodySite = df@sam_data@.Data[[4]])

counts <- samples %>%
  filter(GWColl %in% weeks) %>%
  group_by(SubjectID) %>%
  summarise(N.count = sum(Preg))

SubjectID <- as.character(counts[counts$N.count > 2, ]$SubjectID)
N <- length(SubjectID)

cores <- detectCores()
cl <- makeCluster(cores[1] - 1)
registerDoParallel(cl)
G <- foreach(i = 1:N, .combine = c
             , .packages = c("SpiecEasi", "phyloseq")
) %dopar% {
  set.seed(575)
  otu.subject <- otu[which(sample_data(df)[["SubjectID"]] == SubjectID[i])
                     , -taxa_rare]
  corr <- sparcc(otu.subject)$Cor
  diag(corr) <- 0
  list(corr)
} # end parallel loop
stopCluster(cl)

# Classification ----
Y <- ifelse(Outcomes[Outcomes$SubjectID %in% SubjectID
                     , "Outcome"] == "Term", -1, 1)

source("/R/gp_functions.gibbs_slice.R")
D <- kern.dis((list(G)))
a <- c(0, 2) # flat prior on sigma
b <- c(0, 1)

# Pre-allocate
ns <- 10000
nchains <- 1
N.train <- N - 1 # LOOCV
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

cores <- detectCores()
cl <- makeCluster(cores[1] - 1)
registerDoParallel(cl)
y.prob <- foreach(i = 1:N, .combine = c, .export = "slice"
                  , .packages = "invgamma"
) %dopar% {
  # train/test split
  split_index <- rep("train", N)
  split_index[i] <- "validate"
  
  set.seed(575)
  sims <- mcmc_array(ns, nchains = nchains, params)
  for (ic in 1:nchains)
    sims[, ic, ] <- t(gp.class(dist = D, Y = Y, split_index = split_index
                               , a = a, b = b, ns = ns, monitor = FALSE))
  
  # mixing
  # plot(sims[,1,][, N.train + 1], type = "l")    # sigma
  # plot(sims[,1,][, N.train + 2], type = "l")    # ell
  # plot(sims[,1,][, 1], type = "l")              # f
  
  # predictions
  return(gp.class_pred(sims[,1,], D, split_index, p
                       , burn = .2, thin = 10, avg = TRUE))
} # end for loop
stopCluster(cl)

roc(Y, y.prob, plot = TRUE)
y.pred <- ifelse(y.prob > .5, 1, ifelse(y.prob < .5, -1, NA))
table(Y, y.pred)
mcc(Y, y.pred)
f1(Y, y.pred)

# OCC ----
source("/R/gp_functions.occ.R")

a <- c(1, 1)
b <- c(5, 5)

set.seed(575)
split_index <- ifelse(Y == -1, "train", "validate")
ind.train <- which(split_index == "train")
ind.validate <- which(split_index == "validate")
split_index[sample(ind.validate, 3)] <- "train"
split_index[sample(ind.train, 5)] <- "validate"

ns <- 10000
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

occ.scores <- gp.occ(sims[,1,], burn = .2, thin = 10)

# Figure 8 ----
plot.occ(occ.scores)

# Survival ----
source("/R/gp_functions.survival.R")

X <- as.numeric(as.character(History[History$SubjectID %in% SubjectID, "History"]))
Y <- Outcomes[Outcomes$SubjectID %in% SubjectID, "GDDel"] / 100
D <- kern.dis(list(G), cbind(X, Y))
a <- c(1, 0, 1, 10, 10) # flat prior on Omega and sigma
b <- c(0, 0, 5, 5, 5)

set.seed(575)
ns <- 1000
system.time(test <- gp.survival(D, Y, a, b, N = ns))

# plot(test$Omega, typ = "l")            # Omega
# plot(test$hyperparams[1, ], typ = "l") # sigma
# plot(test$hyperparams[2, ], typ = "l") # ell.net
# plot(test$hyperparams[3, ], typ = "l") # ell.X1
# plot(test$hyperparams[4, ], typ = "l") # ell.time
# plot(sapply(test$l, function(l) l[1]), typ = "l") # l(1)
# plot(sapply(1:dim(test$surv)[3]
#             , function(q) test$surv[1, 1, q]), typ = "l") # l(grid(1), 1)


# Figure 9 ----
plot.surv(test, Y, groups = X, km = T)

# add to plot.surv

# switch names in p1
# geom_line(alpha = .2)
# scale_colour_manual(name = "")

# p1 <- p1 +
#   scale_x_continuous(breaks = c(0, 10, 20, 30, 40) * 7 / 100
#                      , labels = c(0, 10, 20, 30, 40)) +
#   xlab("Gestational Weeks")