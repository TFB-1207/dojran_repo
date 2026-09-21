## Packages
library(dplyr)
library(tidyr)
library(vegan)
library(ggplot2)
library(ggrepel)
library(viridis)

## Load data ##

setwd("C:/Users/kleij057/OneDrive - Wageningen University & Research/DATA/sedaDNA/sedaDNAjune")
plotdna <- read.csv2("cleandna_lowcutoff.csv", sep = ",")

# check data classes
plotdna$depth <- as.numeric(plotdna$depth)
plotdna$rel_abundance <- as.numeric(plotdna$rel_abundance)

## --- 1) Order genus by taxonomy: Phylum -> Order -> Family -> Genus ---

genus_levels <- plotdna %>%
  distinct(order_name, family_name, genus_name) %>%
  arrange(order_name, family_name, genus_name) %>%
  pull(genus_name)

plotdna <- plotdna %>%
  mutate(genus_name = factor(genus_name, levels = genus_levels, ordered = TRUE))


## --- 2) Family ordering and depth as a discrete (equal-height) axis ---

# Family facets ordered by Phylum -> Order -> Family
family_levels <- plotdna %>%
  distinct(order_name, family_name) %>%
  arrange(order_name, family_name) %>%
  pull(family_name)

# Depth factor: deepest at bottom, shallowest at top (so "0 at top")
depth_levels_global <- plotdna %>%
  filter(!is.na(depth)) %>%
  distinct(depth) %>%
  arrange(desc(depth)) %>%      # deepest first
  pull(depth)

plotdna2 <- plotdna %>%
  mutate(
    family_name = factor(family_name, levels = family_levels, ordered = TRUE),
    depth_f     = factor(depth, levels = depth_levels_global, ordered = TRUE)
  )


## Quick, ordered overview heatmap, faceted by family, equal tiles ---

p_pretty <- ggplot(plotdna2, aes(x = genus_name, y = depth_f)) +
  
  # Heatmap tiles, equal-sized in both x and y
  geom_tile(aes(fill = bin), width = 0.9, height = 0.9,
            colour = "black", linewidth = 0, linejoin = "mitre") +
  
  # Color scale
  scale_fill_viridis_d(option = "mako", begin = 0.25, end = 1, name = "Bin") +
  
  # Facet: rows = loc, cols = family
  # scales/space FREE so panels shrink/expand based on data
  facet_grid(
    rows   = vars(loc),
    cols   = vars(family_name),
    scales = "free",     # free x AND y
    space  = "free",     # panel widths & heights adapt to number of taxa/depths
    switch = "y"
  ) +
  
  labs(
    x = "Genus (families ordered by Phylum → Order → Family)",
    y = "Depth",
    title = "Binned abundance heatmap by location and depth, 100 reads min"
  ) +
  
  theme_bw(base_size = 12) +
  theme(
    # Background
    panel.background  = element_rect(fill = "black", colour = NA),
    
    # Grid
    panel.grid.major.x = element_line(colour = "black", linewidth = 0.2),
    panel.grid.minor   = element_line(colour = "black", linewidth = 0.2),
    panel.grid.major.y = element_line(colour = "black", linewidth = 0.2),
    
    # Axes text
    axis.text.x      = element_text(angle = 90, size = 8, vjust = 0.5, hjust = 1, colour = "black"),
    axis.text.y      = element_text(colour = "black"),
    axis.title       = element_text(colour = "black"),
    plot.title       = element_text(face = "bold", colour = "black"),
    
    # Facet strips
    panel.border       = element_rect(colour = "black", linewidth = 0.7, fill = NA),
    strip.background.x = element_rect(fill = "white", colour = "black", linewidth = 0.7),
    strip.background.y = element_rect(fill = "white", colour = "black", linewidth = 0.7),
    strip.text.x       = element_text(size = 7, colour = "black", angle = 90),
    strip.text.y.left  = element_text(angle = 0, face = "bold", colour = "black"),
    strip.placement    = "outside",
    
    panel.spacing.x    = grid::unit(2, "pt"),
    panel.spacing.y    = grid::unit(10, "pt"),
    
    # Legend
    legend.background = element_rect(fill = "white", colour = NA),
    legend.key        = element_rect(fill = "white", colour = NA),
    legend.text       = element_text(colour = "black"),
    legend.title      = element_text(colour = "black")
  )

p_pretty

genera100 <- unique(plotdna2$genus_name)

#################################################################################
# Now run a DCA to inspect the behavior of your samples
# After this we use clustering to identify groups of vegetation that "behave similar"

# remove 116 from clustering becauses the strong persicaria signal hides the rest
plotdna2 <- subset(plotdna2, depth != 116)

# Keep only needed columns
DCA_dat <- plotdna2 %>%
  dplyr::select(loc, depth, genus_name, rel_abundance)

make_comm_matrix <- function(dat,
                             loc_filter = "MO",
                             taxon_col  = "genus_name",
                             depth_col  = "depth",
                             loc_col    = "loc",
                             value_col  = "rel_abundance") {
  
  tax_sym   <- rlang::sym(taxon_col)
  depth_sym <- rlang::sym(depth_col)
  loc_sym   <- rlang::sym(loc_col)
  value_sym <- rlang::sym(value_col)
  
  dat2 <- dat
  
  # Optionally restrict to one location
  if (!is.null(loc_filter)) {
    dat2 <- dat2 %>% dplyr::filter(!!loc_sym == loc_filter)
  }
  
  # Depth must be numeric
  if (!is.numeric(dat2[[rlang::as_string(depth_sym)]])) {
    stop("Depth column is not numeric; convert it before running DCA.")
  }
  
  # Sample metadata
  meta <- dat2 %>%
    dplyr::distinct(!!loc_sym, !!depth_sym) %>%
    dplyr::arrange(!!loc_sym, !!depth_sym) %>%
    dplyr::mutate(sample_id = paste(!!loc_sym, !!depth_sym, sep = "__"))
  
  # Long community data
  comm_long <- dat2 %>%
    dplyr::group_by(!!loc_sym, !!depth_sym, !!tax_sym) %>%
    dplyr::summarise(value = sum(!!value_sym, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(sample_id = paste(!!loc_sym, !!depth_sym, sep = "__"))
  
  # Wide community matrix
  wide <- comm_long %>%
    dplyr::select(sample_id, !!tax_sym, value) %>%
    tidyr::pivot_wider(
      names_from  = !!tax_sym,
      values_from = value,
      values_fill = 0
    )
  
  wide <- wide[match(meta$sample_id, wide$sample_id), ]
  
  if (!identical(wide$sample_id, meta$sample_id)) {
    stop("Row order mismatch between community matrix and metadata.")
  }
  
  cm <- wide %>%
    dplyr::select(-sample_id) %>%
    as.data.frame()
  
  rownames(cm) <- meta$sample_id
  
  # Drop all-zero taxa
  cm <- cm[, colSums(cm, na.rm = TRUE) > 0, drop = FALSE]
  
  list(
    cm   = cm,
    meta = meta
  )
}

run_dca <- function(cm_list) {
  cm   <- cm_list$cm
  meta <- cm_list$meta
  
  dca_obj <- vegan::decorana(cm)
  
  # Site scores
  site_scores <- as.data.frame(vegan::scores(dca_obj, display = "sites"))
  site_scores$sample_id <- rownames(site_scores)
  site_scores <- site_scores %>%
    dplyr::left_join(meta, by = "sample_id")
  
  # Species scores
  species_scores <- as.data.frame(vegan::scores(dca_obj, display = "species"))
  species_scores$taxon <- rownames(species_scores)
  
  list(
    dca     = dca_obj,
    sites   = site_scores,
    species = species_scores
  )
}

# Run DCA
cm_MO  <- make_comm_matrix(DCA_dat, loc_filter = "MO")
dca_MO <- run_dca(cm_MO)

# Sanity check
stopifnot(identical(rownames(cm_MO$cm), dca_MO$sites$sample_id))

# Simple clean DCA plot
p_dca_biplot_MO <- ggplot() +
  geom_point(
    data = dca_MO$sites,
    aes(x = DCA1, y = DCA2, color = depth),
    size = 3
  ) +
  ggrepel::geom_text_repel(
    data = dca_MO$species,
    aes(x = DCA1, y = DCA2, label = taxon),
    color = "black",
    size = 2.8,
    max.overlaps = Inf,
    box.padding = 0.2,
    point.padding = 0.1,
    segment.alpha = 0.4,
    segment.size = 0.2
  ) +
  scale_color_viridis_c(
    option = "magma",
    direction = -1,
    name = "Depth"
  ) +
  labs(
    x = "DCA1",
    y = "DCA2",
    title = "DCA biplot - Madhuvan"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

p_dca_biplot_MO


### Now apply clustering
####################################
# -----------------------------
# 1) Extract DCA species scores
# -----------------------------
library(vegan)
library(factoextra)
taxa_scores <- scores(dca_MO, display = "species")
taxa_scores <- as.data.frame(taxa_scores)

dca_coords <- taxa_scores[, c("DCA1", "DCA2")]

# -----------------------------
# 2) Inspect possible K values
# -----------------------------
fviz_nbclust(dca_coords, kmeans, method = "wss")

# -----------------------------
# 3) Run k-means clustering
# -----------------------------
set.seed(123)

# determine which amount of clusters you want with the broken stick
k <- 4
taxa_clusters <- kmeans(dca_coords, centers = k, nstart = 100)

taxa_scores$cluster <- taxa_clusters$cluster

ggplot(taxa_scores, aes(DCA1, DCA2, color = factor(cluster))) +
  geom_point(size = 3) +
  theme_minimal()

ggplot(taxa_scores, aes(DCA3, DCA4, color = factor(cluster))) +
  geom_point(size = 3) +
  theme_minimal()



# -----------------------------
# 4) Join clusters onto elevation table
#    Assumes 'elevation' has genus_name and q1
#    and that taxon names match genus_name
# -----------------------------
# Ensure cluster is a factor
taxa_scores$cluster <- factor(taxa_scores$cluster)

# Define cluster order (numeric order)
cluster_order <- levels(taxa_scores$cluster)

# Fixed colors for clusters
library(scico)

cluster_cols <- setNames(
  scico(n = length(cluster_order), palette = "managua"),
  cluster_order
)

# -----------------------------
# 6) Create fixed cluster colors
#    Same colors can be reused in all figures
# -----------------------------

cluster_cols

# -----------------------------
# 7) Join cluster onto species labels and fix factor order
# -----------------------------

# Join to species scores (taxon = genus)
species_annot <- dca_MO$species 

species_annot2 <- species_annot %>%
  left_join(
    taxa_scores %>% select(taxon, cluster),
    by = "taxon"
  ) %>%
  mutate(
    cluster = factor(as.character(cluster), levels = cluster_order, ordered = TRUE)
  )

# Also make sure taxa_scores uses the same cluster order
taxa_scores <- taxa_scores %>%
  mutate(
    cluster = factor(as.character(cluster), levels = cluster_order, ordered = TRUE)
  )

# -----------------------------
# 8) Full DCA plot
# -----------------------------

p_dca_biplot_MO <- ggplot() +
  
  # --- Sites (fill = depth) ---
  geom_point(
    data = dca_MO$sites,
    aes(x = DCA1, y = DCA2, fill = depth, shape = loc),
    size = 3.5, colour = "black", alpha = 0.9
  ) +
  scale_shape_manual(
    values = c("MO" = 21, "ST" = 22),
    name = "Location"
  ) +
  scale_fill_viridis(
    option = "mako",
    direction = -1,
    name = "Depth"
  ) +
  
  # Start new fill scale for clusters
  ggnewscale::new_scale_fill() +
  
  # --- Cluster ellipses ---
  stat_ellipse(
    data = species_annot2,
    aes(x = DCA1, y = DCA2, fill = cluster, group = cluster),
    geom = "polygon",
    level = 0.8,
    alpha = 0.25,
    colour = NA
  ) +
  scale_fill_manual(
    values = cluster_cols,
    name = "Cluster",
    drop = FALSE
  ) +
  
  # --- Species labels ---
  ggrepel::geom_text_repel(
    data = species_annot2,
    aes(x = DCA1, y = DCA2, label = taxon, colour = cluster),
    inherit.aes = FALSE,
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.2,
    point.padding = 0.1,
    segment.alpha = 0.4,
    segment.size = 0.2
  ) +
  scale_colour_manual(
    values = cluster_cols,
    name = "Cluster",
    drop = FALSE
  ) +
  
  labs(
    x = "DCA1 51%",
    y = "DCA2 15%",
    title = "DCA biplot - Madhuvan, excluding aquatic plants"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

p_dca_biplot_MO



#######
# -----------------------------
# 8) Full DCA plot
# -----------------------------

p_dca_biplot_MO_34 <- ggplot() +
  
  # --- Sites (fill = depth) ---
  geom_point(
    data = dca_MO$sites,
    aes(x = DCA3, y = DCA4, fill = depth, shape = loc),
    size = 3.5, colour = "black", alpha = 0.9
  ) +
  scale_shape_manual(
    values = c("MO" = 21, "ST" = 22),
    name = "Location"
  ) +
  scale_fill_viridis(
    option = "mako",
    direction = -1,
    name = "Depth"
  ) +
  
  # Start new fill scale for clusters
  ggnewscale::new_scale_fill() +
  
  # --- Cluster ellipses ---
  stat_ellipse(
    data = species_annot2,
    aes(x = DCA3, y = DCA4, fill = cluster, group = cluster),
    geom = "polygon",
    level = 0.8,
    alpha = 0.25,
    colour = NA
  ) +
  scale_fill_manual(
    values = cluster_cols,
    name = "Cluster",
    drop = FALSE
  ) +
  
  # --- Species labels ---
  ggrepel::geom_text_repel(
    data = species_annot2,
    aes(x = DCA3, y = DCA4, label = taxon, colour = cluster),
    inherit.aes = FALSE,
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.2,
    point.padding = 0.1,
    segment.alpha = 0.4,
    segment.size = 0.2
  ) +
  scale_colour_manual(
    values = cluster_cols,
    name = "Cluster",
    drop = FALSE
  ) +
  
  labs(
    x = "DCA3 11%",
    y = "DCA4 8% ",
    title = "DCA biplot - Madhuvan, excluding aquatic plants"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

p_dca_biplot_MO_34

# clusters df
clusterspec <- taxa_scores[,5:6]
write.csv(clusterspec, "clusterslowcutoffk4.csv", row.names = FALSE)

## -----------------------------
## Cluster / DCA inspection block
## -----------------------------

library(vegan)
library(cluster)

## 1) DCA object summary
dca_MO$dca

## 2) Eigenvalues and axis contribution
dca_eig <- eigenvals(dca_MO$dca)

dca_axis_summary <- data.frame(
  axis = names(dca_eig),
  eigenvalue = as.numeric(dca_eig),
  prop = as.numeric(dca_eig / sum(dca_eig)),
  cumulative_prop = as.numeric(cumsum(dca_eig / sum(dca_eig)))
)

dca_axis_summary

## 3) Silhouette based on DCA coordinates
sil_dca <- silhouette(taxa_clusters$cluster, dist(dca_coords))

mean_sil_dca <- mean(sil_dca[, 3])
mean_sil_dca

plot(sil_dca, main = "Silhouette based on DCA coordinates")

## 4) Cluster fit in DCA space
cluster_fit <- taxa_clusters$betweenss / taxa_clusters$totss
cluster_fit

## 5) Silhouette based on original Bray-Curtis taxon distances
taxa_dist <- vegdist(t(cm_MO$cm), method = "bray")

sil_bray <- silhouette(taxa_clusters$cluster, taxa_dist)

mean_sil_bray <- mean(sil_bray[, 3])
mean_sil_bray

plot(sil_bray, main = "Silhouette based on Bray-Curtis taxon distances")

## 6) Mantel test: do DCA distances preserve Bray-Curtis distances?
mantel_test <- mantel(
  taxa_dist,
  dist(dca_coords),
  method = "spearman",
  permutations = 999
)

mantel_test
