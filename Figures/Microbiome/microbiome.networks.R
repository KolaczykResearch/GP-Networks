# Libraries ----
library(SpiecEasi)
library(ggplot2)
library(grid)
library(dplyr)
library(foreach)
library(doParallel)
library(igraph)

# Data ----
load("/Data/PregnancyClosed15.Rdata")

df <- data.frame(SubjectID = PS@sam_data@.Data[[5]]
                 , GDDel = PS@sam_data@.Data[[19]]
                 , GWDel = PS@sam_data@.Data[[20]]
                 , Outcome = PS@sam_data@.Data[[62]]
                 , History = PS@sam_data@.Data[[12]])
df <- unique(df)
Outcomes <- df[, c("SubjectID", "Outcome", "GWDel", "GDDel")]
History <- df[, c("SubjectID", "History")]

# Histogram ----
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
df <- subset_samples(df, BodySite == "Vaginal_Swab")

otu <- otu_table(df)@.Data
otu.present <- otu; otu.present[otu > 0] <- 1
taxa_rare <- which(colSums(otu.present) < 10 | colSums(otu) < 500)

samples <- data.frame(SubjectID = df@sam_data@.Data[[5]]
                      , GDColl = df@sam_data@.Data[[8]]
                      , GWColl = df@sam_data@.Data[[9]]
                      , Preg = df@sam_data@.Data[[56]]
                      , BodySite = df@sam_data@.Data[[4]])

counts <- samples %>%
  filter(GWColl %in% weeks) %>%
  filter(BodySite == "Vaginal_Swab") %>%
  group_by(SubjectID) %>%
  summarise(m.count = sum(Preg))

SubjectID <- as.character(counts[counts$m.count > 1, ]$SubjectID)
m <- length(SubjectID)

cores <- detectCores()
cl <- makeCluster(cores)
registerDoParallel(cl)
G <- foreach(k = 1:m, .combine = c
             , .packages = c("SpiecEasi", "phyloseq")
) %dopar% {
  set.seed(575)
  otu.subject <- otu[which(sample_data(df)[["SubjectID"]] == SubjectID[k])
                     , -taxa_rare]
  corr <- sparcc(otu.subject)$Cor
  diag(corr) <- 0
  list(corr)
} # end parallel loop
stopCluster(cl)

Y <- ifelse(Outcomes[Outcomes$SubjectID %in% SubjectID
                     , "Outcome"] == "Term", -1, 1)

X <- as.numeric(as.character(History[History$SubjectID %in% SubjectID, "History"]))
GDDel <- Outcomes[Outcomes$SubjectID %in% SubjectID, "GDDel"] / 100

# Visualization ----
taxa <- tax_table(df)@.Data
taxa.genus <- taxa[-taxa_rare, "Genus"]
taxa.species <- taxa[-taxa_rare, "Species"]
taxa.phylum <- taxa[-taxa_rare, "Phylum"]

otu.prop <- apply(otu, 2, function(x) x / sum(x))

par(mfrow = c(1, 2))
par(mar = c(0, 0, 0, 0) + .1)

net1 <- graph_from_adjacency_matrix(G[[19]]
                                    , mode = "undirected", weighted = TRUE)
otu1 <- otu[which(sample_data(df)[["SubjectID"]] == SubjectID[19])
            , -taxa_rare]

V(net1)$taxa <- colMeans(otu1)
V(net1)$color <- as.factor(taxa.phylum)
V(net1)$label <- NA

set.seed(575)
layout <- layout_with_dh(net1)
plot(net1
     , edge.width = E(net1)$weight
     , edge.color = ifelse(E(net1)$weight < 0, "gray30", "orangered3")
     , vertex.size = V(net1)$taxa^(1/3)*2
     , vertex.label.font = 4
     , vertex.label.color = "gray30", vertex.label.cex = 1
     , layout = layout
     , xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2), asp = 0)
title("Preterm", line = -14.5)

net2 <- graph_from_adjacency_matrix(G[[10]]
                                    , mode = "undirected", weighted = TRUE)
otu2 <- otu[which(sample_data(df)[["SubjectID"]] == SubjectID[10])
            , -taxa_rare]

V(net2)$taxa <- colMeans(otu2)
V(net2)$color <- as.factor(taxa.phylum)
V(net2)$label <- NA

plot(net2
     , edge.width = E(net2)$weight
     , edge.color = ifelse(E(net2)$weight < 0, "gray30", "orangered3")
     , vertex.size = V(net2)$taxa^(1/3)*2
     , vertex.label.font = 4
     , vertex.label.color = "gray30", vertex.label.cex = 1
     , layout = layout
     , xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2), asp = 0)
title("Term", line = -14.5)

title("Individual Microbiome Networks", line = -1.5, outer = TRUE)

# EDA ----
# need to use absolute value but doesn't distinguish between
#    positive and negative weights
# could binarize but wouldn't capture strength of edge

G_abs <- lapply(G, abs)
G_inv <- lapply(G_abs, function(x) 1/x)
G_net <- lapply(G_abs, graph_from_adjacency_matrix, mode = "undirected", weighted = TRUE)
G_dis <- lapply(G_inv, graph_from_adjacency_matrix, mode = "undirected", weighted = TRUE)

# Network statistics ----

# Total strength
strength <- sapply(G_net, strength)
t.test(colSums(strength[, Y == -1]), colSums(strength[, Y == 1]))
sd(colSums(strength[, Y == -1])); sd(colSums(strength[, Y == 1]))

# Number of communities
community_num <- sapply(G_net, function(g) length(fastgreedy.community(g)))
t.test(community_num[Y == -1], community_num[Y == 1])
sd(community_num[Y == -1]); sd(community_num[Y == 1])

# Closeness
v_close <- sapply(G_dis, closeness)
t.test(colSums(v_close[, Y == -1]), colSums(v_close[, Y == 1]))
sd(colSums(v_close[, Y == -1])); sd(colSums(v_close[, Y == 1]))

# Diameter
diam <- sapply(G_dis, function(g) diameter(g, directed = FALSE))
t.test(diam[Y == -1], diam[Y == 1])
sd(diam[Y == -1]); sd(diam[Y == 1])

# Node statistics ----
ind.lact <- grep("Lactobacillus", taxa.species, fixed = FALSE)
ind.gard <- grep("Gardnerella", taxa.species, fixed = FALSE)
ind.urea <- grep("Ureaplasma", taxa.species, fixed = FALSE)
ind.fin <- grep("Finegoldia", taxa.species, fixed = FALSE)

# Central nodes by strength
mean(apply(strength[, Y == -1], 2, function(x) ind.gard %in% order(x, decreasing = TRUE)[1:3]))
mean(apply(strength[, Y == 1], 2, function(x) ind.gard %in% order(x, decreasing = TRUE)[1:3]))

mean(apply(strength[, Y == -1], 2, function(x) ind.urea %in% order(x, decreasing = TRUE)[1:3]))
mean(apply(strength[, Y == 1], 2, function(x) ind.urea %in% order(x, decreasing = TRUE)[1:3]))

mean(apply(strength[, Y == -1], 2, function(x) any(order(x, decreasing = TRUE)[1:3] %in% ind.lact)))
mean(apply(strength[, Y == 1], 2, function(x) any(order(x, decreasing = TRUE)[1:3] %in% ind.lact)))

mean(apply(strength[, Y == -1], 2, function(x) ind.fin %in% order(x, decreasing = TRUE)[1:3]))
mean(apply(strength[, Y == 1], 2, function(x) ind.fin %in% order(x, decreasing = TRUE)[1:3]))

# Central nodes by closeness
mean(apply(v_close[, Y == -1], 2, function(x) ind.gard %in% order(x, decreasing = TRUE)[1:3]))
mean(apply(v_close[, Y == 1], 2, function(x) ind.gard %in% order(x, decreasing = TRUE)[1:3]))

mean(apply(v_close[, Y == -1], 2, function(x) ind.urea %in% order(x, decreasing = TRUE)[1:3]))
mean(apply(v_close[, Y == 1], 2, function(x) ind.urea %in% order(x, decreasing = TRUE)[1:3]))

mean(apply(v_close[, Y == -1], 2, function(x) any(order(x, decreasing = TRUE)[1:3] %in% ind.lact)))
mean(apply(v_close[, Y == 1], 2, function(x) any(order(x, decreasing = TRUE)[1:3] %in% ind.lact)))

mean(apply(v_close[, Y == -1], 2, function(x) ind.fin %in% order(x, decreasing = TRUE)[1:3]))
mean(apply(v_close[, Y == 1], 2, function(x) ind.fin %in% order(x, decreasing = TRUE)[1:3]))

# Central nodes by eigencentrality
eigen_cen <- sapply(G_net, function(g) eigen_centrality(g)$vector)

mean(apply(eigen_cen[, Y == -1], 2, function(x) ind.gard %in% order(x, decreasing = TRUE)[1:3]))
mean(apply(eigen_cen[, Y == 1], 2, function(x) ind.gard %in% order(x, decreasing = TRUE)[1:3]))

mean(apply(eigen_cen[, Y == -1], 2, function(x) ind.urea %in% order(x, decreasing = TRUE)[1:3]))
mean(apply(eigen_cen[, Y == 1], 2, function(x) ind.urea %in% order(x, decreasing = TRUE)[1:3]))

mean(apply(eigen_cen[, Y == -1], 2, function(x) any(order(x, decreasing = TRUE)[1:3] %in% ind.lact)))
mean(apply(eigen_cen[, Y == 1], 2, function(x) any(order(x, decreasing = TRUE)[1:3] %in% ind.lact)))

mean(apply(eigen_cen[, Y == -1], 2, function(x) ind.fin %in% order(x, decreasing = TRUE)[1:3]))
mean(apply(eigen_cen[, Y == 1], 2, function(x) ind.fin %in% order(x, decreasing = TRUE)[1:3]))

# Central nodes by hub score
hub_cen <- sapply(G_net, function(g) hub_score(g)$vector)
mean(apply(hub_cen[, Y == -1], 2, function(x) ind.gard %in% order(x, decreasing = TRUE)[1:3]))
mean(apply(hub_cen[, Y == 1], 2, function(x) ind.gard %in% order(x, decreasing = TRUE)[1:3]))

mean(apply(hub_cen[, Y == -1], 2, function(x) ind.urea %in% order(x, decreasing = TRUE)[1:3]))
mean(apply(hub_cen[, Y == 1], 2, function(x) ind.urea %in% order(x, decreasing = TRUE)[1:3]))

mean(apply(hub_cen[, Y == -1], 2, function(x) any(order(x, decreasing = TRUE)[1:3] %in% ind.lact)))
mean(apply(hub_cen[, Y == 1], 2, function(x) any(order(x, decreasing = TRUE)[1:3] %in% ind.lact)))

mean(apply(hub_cen[, Y == -1], 2, function(x) ind.fin %in% order(x, decreasing = TRUE)[1:3]))
mean(apply(hub_cen[, Y == 1], 2, function(x) ind.fin %in% order(x, decreasing = TRUE)[1:3]))
