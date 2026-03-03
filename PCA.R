
library(tidyverse)

overlap_mat <- read.table("Data Files/Orthogroups_SpeciesOverlaps.tsv",
                      header = TRUE,
                      row.names = 1,
                      sep = "\t",
                      check.names = FALSE)

# Read orthogroup totals
og_counts <- read.csv("Data Files/orthogroup_counts.csv")


family_map <- read.csv("Data Files/Trait Data/Families.csv")
diet_map <- read.csv("Data Files/Trait Data/Diet.csv")
nesting_map <- read.csv("Data Files/Trait Data/Nesting_Substrate.csv")
sociality_map <- read.csv("Data Files/Trait Data/Sociality.csv")

# ----------------------------
# NORMALIZATION STEP
# ----------------------------
# Scale overlap values to 0–1 with Jaccard distances

# Make named vector
og_vec <- og_counts$Total_Orthogroups
names(og_vec) <- og_counts$Species

# Ensure species order matches matrix
og_vec <- og_vec[rownames(overlap_mat)]

# Normalization by total number of orthogroups per species
min_mat <- outer(og_vec, og_vec, pmin)
norm_mat2 <- overlap_mat / min_mat

overlap_mat <- norm_mat2

#remove outliers and ref species
overlap_mat <- overlap_mat[!rownames(overlap_mat) %in% c("Drosophila_melanogaster", "Anopheles_gambiae"), ]

# Convert similarity to distance
dist_mat <- 1 - norm_mat2
dist_obj <- as.dist(dist_mat)


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

# Variance explained
var_exp <- round(pcoa$eig / sum(pcoa$eig) * 100, 1)

# ----------------------------
# Plot
# ----------------------------
ggplot(scores2, aes(PCoA1, PCoA2, color = Family, label = Species)) +
  geom_point(size = 2) +
  #geom_text(vjust = -0.6, size = 3) +
  theme_bw() +
  labs(
    title = "PCoA of Gene Overlap by Family",
    x = paste0("PCoA1 (", var_exp[1], "%)"),
    y = paste0("PCoA2 (", var_exp[2], "%)"),
    color = "Diet"
  )

#Zoom in on clusters
ggplot(scores2, aes(PCoA1, PCoA2, color = Family)) +
  geom_point(size = 2) +
  coord_cartesian(xlim = c(-0.03, 0.04), ylim = c(-0.06, 0.04)) +
  theme_bw() +
  labs(
    title = "PCoA of Gene Overlap by Family",
    x = paste0("PCoA1 (", var_exp[1], "%)"),
    y = paste0("PCoA2 (", var_exp[2], "%)"),
    color = "Family"
  )

#This method tells us how similar species are based on the number of genes they share


# ----------------------------
# PERMANOVA
# ----------------------------

library(vegan)
# Distance matrix in ordination space
ord_dist <- dist(scores2[, c("PCoA1", "PCoA2")])

# PERMANOVA
adonis_result <- adonis2(ord_dist ~ Family, data = scores2, permutations = 999)

adonis_result

