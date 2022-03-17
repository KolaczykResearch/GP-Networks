# Version : R 4.0.5
# Libraries ----
library(ggplot2)
library(dplyr)
library(igraph)
library(reshape2)

# Data ----
load("/Data/PregnancyClosed15.Rdata")

weeks <- 1:30
df <- subset_samples(PS, GWColl %in% weeks)
df <- subset_samples(df, BodySite == "Vaginal_Swab")

# subjects
samples <- data.frame(SubjectID = df@sam_data@.Data[[5]]
                      , GDColl = df@sam_data@.Data[[8]]
                      , GWColl = df@sam_data@.Data[[9]]
                      , Preg = df@sam_data@.Data[[56]]
                      , BodySite = df@sam_data@.Data[[4]]
                      , Outcome = df@sam_data@.Data[[62]])

counts <- samples %>%
  group_by(SubjectID) %>%
  summarise(m.count = sum(Preg)
            , Outcome = unique(Outcome))

SubjectID <- as.character(counts[counts$m.count > 1, ]$SubjectID)

# OTUs
otu <- otu_table(df)@.Data

# taxa
taxa <- tax_table(df)@.Data
taxa.species <- taxa[, "Species"]

ind.lact <- grep("Lactobacillus", taxa.species, fixed = FALSE)
ind.gard <- grep("Gardnerella", taxa.species, fixed = FALSE)
ind.urea <- grep("Ureaplasma", taxa.species, fixed = FALSE)
ind.fin <- grep("Finegoldia", taxa.species, fixed = FALSE)

# OCC results
# load("./Data/occ_lambda.RData")

# networks
G.abs <- lapply(G[split_index == "train"], function(g) {
  graph_from_adjacency_matrix(abs(g)
                              , mode = "undirected"
                              , weighted = TRUE)
})
G.dis <- lapply(G[split_index == "train"], function(g) {
  graph_from_adjacency_matrix(1/abs(g)
                              , mode = "undirected"
                              , weighted = TRUE)
})

# Post-hoc ----
f.train <- apply(sims[, 1, seq_len(m.train)], 2, mean)

df_post_hoc <- data.frame(Latent = f.train
                          , Outcome = counts[which(split_index == "train"), "Outcome"]
                          , Lactobacillus = sapply(SubjectID[which(split_index == "train")]
                                                   , function(k) mean(otu[which(sample_data(df)[["SubjectID"]] == k)
                                                                          , ind.lact]))
                          , Gardnerella = sapply(SubjectID[which(split_index == "train")]
                                                 , function(k) mean(otu[which(sample_data(df)[["SubjectID"]] == k)
                                                                        , ind.gard]))
                          , Strength = colSums(sapply(G.abs, function(g) strength(g)))
                          , Closeness = colSums(sapply(G.dis, function(g) closeness(g)))
                          , Diameter = sapply(G.dis, function(g) diameter(g, directed = FALSE))
                          , Betweenness = colSums(sapply(G.dis, function(g) betweenness(g))))
df_post_hoc$Outcome <- factor(df_post_hoc$Outcome
                              , levels = c("VeryPreterm", "Marginal", "Term")
                              , labels = c("Preterm", "Preterm", "Term"))

df_plot <- melt(subset(df_post_hoc, select = c("Latent", "Outcome", "Lactobacillus", "Closeness", "Diameter"))
                    , id.vars = c("Latent", "Outcome"))

ggplot(df_plot, aes(x = Latent, y = value, color = Outcome)) +
  geom_point() +
  theme_bw() +
  theme(panel.spacing = unit(1, "lines")) +
  ylab(NULL) +
  facet_grid(variable~., scales = "free_y")