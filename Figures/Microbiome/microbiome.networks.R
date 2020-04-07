# Libraries ----
library(SpiecEasi)
library(grid)
library(dplyr)
library(foreach)
library(doParallel)
library(igraph)

# Data ----
# from : http://statweb.stanford.edu/~susan/papers/PNASRR.html
load("/R/PregnancyClosed15.Rdata")

df <- data.frame(SubjectID = PS@sam_data@.Data[[5]]
                 , GDDel = PS@sam_data@.Data[[19]]
                 , GWDel = PS@sam_data@.Data[[20]]
                 , Outcome = PS@sam_data@.Data[[62]]
                 , History = PS@sam_data@.Data[[12]])
df <- unique(df)
Outcomes <- df[, c("SubjectID", "Outcome", "GWDel", "GDDel")]
History <- df[, c("SubjectID", "History")]

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

Y <- ifelse(Outcomes[Outcomes$SubjectID %in% SubjectID
                     , "Outcome"] == "Term", -1, 1)

# Taxa descriptors ----
df.vagina <- subset_samples(PS, BodySite == "Vaginal_Swab" & GWColl %in% weeks)
otu.vagina <- otu_table(df.vagina)@.Data
otu.vagina <- otu.vagina[, -taxa_rare]

df.saliva <- subset_samples(PS, BodySite == "Saliva" & GWColl %in% weeks)
otu.saliva <- otu_table(df.saliva)@.Data
otu.saliva <- otu.saliva[, -taxa_rare]

df.stool <- subset_samples(PS, BodySite == "Stool" & GWColl %in% weeks)
otu.stool <- otu_table(df.stool)@.Data
otu.stool <- otu.stool[, -taxa_rare]

df.tooth <- subset_samples(PS, BodySite == "Tooth_Gum" & GWColl %in% weeks)
otu.tooth <- otu_table(df.tooth)@.Data
otu.tooth <- otu.tooth[, -taxa_rare]

otu.prev <- cbind(colSums(otu.vagina), colSums(otu.saliva)
                  , colSums(otu.stool), colSums(otu.tooth))
otu.prop <- apply(otu.prev, 2, function(x) x / sum(x))
site <- apply(otu.prop, 1, which.max)

# Figure 2 ----
par(mfrow = c(1, 2))
par(mar = c(0, 0, 0, 0) + .1)

net1 <- graph_from_adjacency_matrix(G[[10]]
                                    , mode = "undirected", weighted = TRUE)
otu1 <- otu[which(sample_data(df)[["SubjectID"]] == SubjectID[10])
            , -taxa_rare]

V(net1)$taxa <- colMeans(otu1)
V(net1)$site <- ifelse(site == 1, "lightpink", ifelse(site == 2, "dodgerblue4"
                                                      , ifelse(site == 3, "goldenrod2", "chartreuse4")))
V(net1)$label <- NA

set.seed(575)
layout <- layout_with_dh(net1)
plot(net1
     , edge.width = E(net1)$weight^2
     , edge.color = ifelse(E(net1)$weight < 0, "gray30", "orangered3")
     , vertex.size = V(net1)$taxa^(1/3)*2
     , vertex.color = V(net1)$site
     , vertex.label.font = 4
     , vertex.label.color = "gray30", vertex.label.cex = 1
     , layout = layout
     , xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2), asp = 0)
title("Preterm", line = -14.5)

net2 <- graph_from_adjacency_matrix(G[[18]]
                                    , mode = "undirected", weighted = TRUE)
otu2 <- otu[which(sample_data(df)[["SubjectID"]] == SubjectID[18])
            , -taxa_rare]

V(net2)$taxa <- colMeans(otu2)
V(net2)$site <- ifelse(site == 1, "lightpink", ifelse(site == 2, "dodgerblue4"
                                                      , ifelse(site == 3, "goldenrod2", "chartreuse4")))
V(net2)$label <- NA

plot(net2
     , edge.width = E(net2)$weight^2
     , edge.color = ifelse(E(net2)$weight < 0, "gray30", "orangered3")
     , vertex.size = V(net2)$taxa^(1/3)*2
     , vertex.color = V(net2)$site
     , vertex.label.font = 4
     , vertex.label.color = "gray30", vertex.label.cex = 1
     , layout = layout
     , xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2), asp = 0)
title("Term", line = -14.5)

title("Individual Microbiome Networks", line = -1.5, outer = TRUE)