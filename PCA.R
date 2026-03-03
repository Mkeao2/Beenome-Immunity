library(tidyverse)

overlap <- read.table("Orthogroups_SpeciesOverlaps.tsv",
                      header = TRUE,
                      row.names = 1,
                      sep = "\t",
                      check.names = FALSE)

family_map <- read.csv("species_family - Sheet1.csv")
diet_map <- read.csv("Trait data  - Diet.csv")
nesting_map <- read.csv("Trait data  - Nesting_substrate (1).csv")
sociality_map <- read.csv("Trait data  - Sociality (1).csv")

# Convert to numeric matrix
overlap_mat <- as.matrix(overlap)
overlap_mat <- apply(overlap_mat, 2, as.numeric)
rownames(overlap_mat) <- rownames(overlap)

# Optional sanity check
str(overlap_mat)

#remove outliers and ref species
overlap_mat <- overlap_mat[!rownames(overlap_mat) %in% c("Drosophila_melanogaster", "Anopheles_gambiae"), ]

# ----------------------------
# NORMALIZATION STEP
# ----------------------------
# Scale overlap values to 0–1
max_overlap <- max(overlap_mat, na.rm = TRUE)
overlap_norm <- overlap_mat / max_overlap

# Convert similarity to distance
dist_mat <- 1 - overlap_norm
dist_obj <- as.dist(dist_mat)


# ----------------------------
# PCoA
# ----------------------------
pcoa <- cmdscale(dist_obj, eig = TRUE, k = 2)

# Extract scores
scores <- as.data.frame(pcoa$points)
scores$Species <- rownames(scores)
colnames(scores)[1:2] <- c("PCoA1", "PCoA2")

scores2 <- left_join(scores, diet_map, by = "Species")

string_to_remove <- "NA"
scores2 <- scores2 %>%
  filter(!grepl(string_to_remove, Diet_broad))

# Variance explained
var_exp <- round(pcoa$eig / sum(pcoa$eig) * 100, 1)

# ----------------------------
# Plot
# ----------------------------
ggplot(scores2, aes(PCoA1, PCoA2, color = Diet_broad, label = Species)) +
  geom_point(size = 2) +
  #geom_text(vjust = -0.6, size = 3) +
  theme_bw() +
  stat_ellipse(level = 0.95) +
  labs(
    title = "PCoA of Gene Overlap by Diet",
    x = paste0("PCoA1 (", var_exp[1], "%)"),
    y = paste0("PCoA2 (", var_exp[2], "%)"),
    color = "Diet"
  )

#Zoom in on clusters
ggplot(scores2, aes(PCoA1, PCoA2, color = Family)) +
  geom_point(size = 2) +
  coord_cartesian(xlim = c(-0.06, 0.1), ylim = c(-0.1, 0.1)) +
  theme_bw() +
  stat_ellipse(level = 0.95) +
  labs(title = "Zoom on most similar species region")

#This method tells us how similar species are based on the number of genes they share


# ----------------------------
# PERMANOVA
# ----------------------------

# Distance matrix in ordination space
ord_dist <- dist(scores2[, c("PCoA1", "PCoA2")])

# PERMANOVA
adonis_result <- adonis2(ord_dist ~ Diet_broad, data = scores2, permutations = 999)

adonis_result

