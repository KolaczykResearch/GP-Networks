# Version : R 4.0.5
# Libraries ----
library(vegan)
library(lme4)
library(dplyr)
library(kernInt)

# Functions ----
source(".R/gp_functions.gibbs_slice.R")

# Data ----
load(".Data/PregnancyClosed15.Rdata")

weeks <- 1:30
df <- subset_samples(PS, GWColl %in% weeks)
df <- subset_samples(df, BodySite == "Vaginal_Swab")

otu <- otu_table(df)@.Data

samples <- data.frame(SubjectID = df@sam_data@.Data[[5]]
                      , GDColl = df@sam_data@.Data[[8]]
                      , GWColl = df@sam_data@.Data[[9]]
                      , Preg = df@sam_data@.Data[[56]]
                      , BodySite = df@sam_data@.Data[[4]]
                      , Outcome = df@sam_data@.Data[[62]]
                      , History = df@sam_data@.Data[[12]])

subjects <- unique(samples$SubjectID)
m <- length(subjects)

df.comp <- samples[, c("SubjectID", "GDColl")]
df.comp$preterm <- ifelse(samples$Outcome == "Term", 0, 1)

# Diversity ----
df.comp$alpha <- apply(otu, 1, vegan::diversity, index = "shannon")
for (k in seq_len(m)) {
  df.comp[df.comp$SubjectID == subjects[k]
          , "beta"] <- mean(vegdist(otu[which(sample_data(df)[["SubjectID"]] == subjects[k]), ],
                                    binary = TRUE))
}

# kernInt ----
df.kern <- otu
rownames(df.kern) <- df.comp$SubjectID

# LOOCV ----
set.seed(575)

results.alpha <- results.beta <- results.kernInt <- vector(length = m)

for (k in seq_len(m)) {
  
  data_train <- df.comp[df.comp$SubjectID != subjects[k], ]
  data_test <- df.comp[df.comp$SubjectID == subjects[k], ]
  
  # https://lme4.r-forge.r-project.org/book/Ch4.pdf
  mm <- glmer(preterm ~ 1 + alpha*GDColl + (1 | SubjectID) + (0 + GDColl | SubjectID)
              , family = binomial
              , data = data_train)
  results.alpha[k] <- mean(predict(mm, newdata = data_test
                                   , allow.new.levels = TRUE
                                   , type = "response"))
  
  mm <- glm(preterm ~ beta
            , family = "binomial"
            , data = unique(data_train[, c("preterm", "beta")]))
  results.beta[k] <- predict(mm, newdata = unique(data_test[, c("preterm", "beta")])
                             , type = "response")
  
  
  mm <- classify(data = df.kern
                 , y = df.comp$preterm
                 , kernel = "clin"
                 , p = which(rownames(df.kern) == subjects[k])
                 , prob = TRUE)
  results.kernInt[k] <- mean(mm$prediction[, "1"])
  
}

# Results ----
outcomes <- df.comp %>% group_by(SubjectID) %>% summarise(Y = mean(preterm))
Y <- outcomes$Y

sum(ifelse(results.alpha > .5, 1, 0) == Y) / m
roc(labels = Y, scores = results.alpha)
mcc(Y * 2 - 1, ifelse(results.alpha > .5, 1, -1))
f1(Y * 2 - 1, ifelse(results.alpha > .5, 1, -1))

sum(ifelse(results.beta > .5, 1, 0) == Y) / m
roc(labels = Y, scores = results.beta)
mcc(Y * 2 - 1, results.beta)
f1(Y * 2 - 1, results.beta)

sum(ifelse(results.kernInt > .5, 1, 0) == Y) / m
roc(labels = Y, scores = results.kernInt)
mcc(Y * 2 - 1, results.kernInt)
f1(Y * 2 - 1, results.kernInt)
