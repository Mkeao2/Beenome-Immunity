
library(tidyverse)
library(dplyr)
library(readr)
library(vegan)


# Read orthogroup totals
og_counts <- read.csv("Data Files/orthogroup_counts.csv")
orthogroups <- read_tsv("Data Files/Orthogroups.tsv")

family_map <- read.csv("Data Files/Trait Data/Families.csv")
diet_map <- read.csv("Data Files/Trait Data/Diet.csv")
nesting_map <- read.csv("Data Files/Trait Data/Nesting_Substrate.csv")
sociality_map <- read.csv("Data Files/Trait Data/Sociality_Sarah.csv")

# -------------------------------------------
# Option to filter to only include immune orthogroups
# -------------------------------------------

immune_og <- read_csv("Data Files/Immune_Orthogroup_Annotation_Current.csv")

immune_ids <- immune_og$Orthogroup

# FOR IMMUNE ANALYSIS ONLY! Filter Orthogroups file
immune_orthogroups <- orthogroups %>%
  filter(Orthogroup %in% immune_ids)

# -------------------------------------------
# Transform counts to presence/absence matrix
# -------------------------------------------

presence_mat <- orthogroups

presence_mat[,-1] <- ifelse(is.na(presence_mat[,-1]), 0, 1)

presence_mat <- as.matrix(presence_mat[,-1])

#remove outliers and ref species
overlap_mat <- overlap_mat[!rownames(overlap_mat) %in% c("Drosophila_melanogaster", "Anopheles_gambiae"), ]

#Jaccard distance
dist_obj <- vegdist(t(presence_mat), method = "jaccard")

t(presence_mat)


# ----------------------------
# PCoA
# ----------------------------
pcoa <- cmdscale(dist_obj, eig = TRUE, k = 2)

# Extract scores
scores <- as.data.frame(pcoa$points)
scores$Species <- rownames(scores)
colnames(scores)[1:2] <- c("PCoA1", "PCoA2")

scores2 <- left_join(scores, family_map, by = "Species")

#removes rows with missing values
scores2 <- scores2 %>% filter(!is.na(Family))

string_to_remove <- "Semisocial" 
scores2 <- scores2 %>% filter(!grepl(string_to_remove, Sociality))

# Variance explained
var_exp <- round(pcoa$eig / sum(pcoa$eig) * 100, 1)

# ----------------------------
# PERMANOVA
# ----------------------------
# Distance matrix in ordination space
ord_dist <- dist(scores2[, c("PCoA1", "PCoA2")])

# PERMANOVA
adonis2(ord_dist ~ Family, data = scores2, permutations = 999)

# ----------------------------
# Plot
# ----------------------------
  
ggplot(scores2, aes(PCoA1, PCoA2, color = Family, label = Family)) +
  geom_point(size = 2) +
  #geom_text(vjust = -0.6, size = 3) +
  theme_bw() +
  labs(
    title = "PCoA of Gene Overlap by Family - All Orthogroups",
    x = paste0("PCoA1 (", var_exp[1], "%)"),
    y = paste0("PCoA2 (", var_exp[2], "%)"),
    color = "Sociality"
  ) + 
  annotate("text", x = 0.05, y = 0.05,
           label = "R^2 == 0.02", 
           parse = TRUE) +
  annotate("text", x = 0.05, y = 0.037,
           label = "P == 0.56", 
           parse = TRUE)
           hjust = 0.1, vjust = 1)


#Zoom in on clusters
ggplot(scores2, aes(PCoA1, PCoA2, color = Family)) +
  geom_point(size = 2) +
  coord_cartesian(xlim = c(-0.02, 0.01), ylim = c(-0.05, 0.04)) +
  theme_bw() +
  labs(
    title = "PCoA of Gene Overlap by Sociality",
    x = paste0("PCoA1 (", var_exp[1], "%)"),
    y = paste0("PCoA2 (", var_exp[2], "%)"),
    color = "Sociality"
  ) 




