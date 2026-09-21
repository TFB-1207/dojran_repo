# steps

library(ggplot2)
#install.packages("ggridges")
library(ggridges)

library(dplyr)
library(tidyr)
library(plotly)

# Load data
#setwd("C:/projects/repositories/zo_repo/260303_taxa_out.tsv")
nepaldat <- readr::read_tsv("C:/projects/repositories/zo_repo/260303_taxa_out.tsv")
nepaldat <- nepaldat[,c(1:4,7:70)]
labels <- read.csv2("C:/projects/repositories/zo_repo/Lib_id2sample_id_new.csv", sep = ",")

# Rename the column
names(nepaldat)[names(nepaldat) == "library"] <- "Library_id"

# Reformat the IDs
nepaldat$Library_id <- gsub("_", ".", nepaldat$Library_id)

#add informative ID instead of library ID
#colnames(labels)[1] <- "lib_id"
relab <- merge(labels, nepaldat, by = "Library_id")

# Now we remove negative controls
dat_exneg <- relab[!relab$Sample_id == "Extraction negative control",]
dat_exneg <- dat_exneg[!dat_exneg$Sample_id == "Library negative control",]

labels2 <- dat_exneg$Sample_id
labels2 <- as.data.frame(labels2)

labels2$loc <- sub("[0-9]+", "", labels2$labels2)
labels2$depth <- as.numeric(sub("[A-Za-z]+", "", labels2$labels2))

dat_exneg$Sample_id <- labels2$depth
names(dat_exneg)[names(dat_exneg) == "Sample_id"] <- "depth"
names(dat_exneg)[names(dat_exneg) == "Site_name"] <- "loc"
dat_exneg$loc <- labels2$loc

sedadna <- dat_exneg



# Create plotting DF
# filter on relevant phyla
adna <- subset(sedadna,  phylum_name == "Streptophyta")

#only more than X readcounts
adna <- adna[adna$genus_count_assembly > 200,]

# ony keep genera that occur in >3 samples
adna_subset <- adna %>%
  group_by(genus_name) %>%
  filter(n_distinct(depth) >= 3)

adna <- adna_subset

#only if occurs in 3 or more samples

# make label per observation
# taxa+site+depth

adna$identifier<- paste(adna$depth, adna$loc, adna$taxa_name, sep = "_")

# now subset to simple df

adna <- adna[,c(71,2,3,17,19,38,41)]


adnaMO <- subset(adna, adna$loc =="MO")
adnaST<- subset(adna, adna$loc =="ST")


# 1. Reshape to long format: one row per depth × measure
plot_dat_MO <- adnaMO %>%
  select(depth, identifier, deam3p_frac_genus, deam5p_frac_genus) %>%
  pivot_longer(
    cols = c(deam3p_frac_genus, deam5p_frac_genus),
    names_to = "end",
    values_to = "frac"
  )

# 2. Make interactive scatter plot
figMO <- plot_ly(
  data = plot_dat_MO,
  x = ~depth,
  y = ~frac,
  type = "scatter",
  mode = "markers",
  color = ~end,  # one color for each column
  text = ~paste(
    "identifier:", identifier,
    "<br>depth:", depth,
    "<br>", end, ":", round(frac, 3)
  ),
  hoverinfo = "text"
) %>%
  plotly::layout(
    xaxis = list(title = "Depth"),
    yaxis = list(title = "Deamination fraction (genus)"),
    legend = list(title = list(text = "End"))
  )

figMO


# Long format
plot_dat_MO <- adnaMO %>%
  select(depth, identifier, deam3p_frac_genus, deam5p_frac_genus) %>%
  pivot_longer(
    cols = c(deam3p_frac_genus, deam5p_frac_genus),
    names_to = "end",
    values_to = "frac"
  )


# 1. Reshape to long format: one row per depth × measure
plot_dat_ST <- adnaST %>%
  select(depth, identifier, deam3p_frac_genus, deam5p_frac_genus) %>%
  pivot_longer(
    cols = c(deam3p_frac_genus, deam5p_frac_genus),
    names_to = "end",
    values_to = "frac"
  )

# 2. Make interactive scatter plot
figST <- plot_ly(
  data = plot_dat_ST,
  x = ~depth,
  y = ~frac,
  type = "scatter",
  mode = "markers",
  color = ~end,  # one color for each column
  text = ~paste(
    "identifier:", identifier,
    "<br>depth:", depth,
    "<br>", end, ":", round(frac, 3)
  ),
  hoverinfo = "text"
) %>%
  layout(
    xaxis = list(title = "Depth"),
    yaxis = list(title = "Deamination fraction (genus)"),
    legend = list(title = list(text = "End"))
  )

figST


### weighted ###
plot_dat_MO <- adnaMO %>%
  select(depth, identifier, genus_count_assembly,
         deam3p_frac_genus, deam5p_frac_genus) %>%
  pivot_longer(
    cols = c(deam3p_frac_genus, deam5p_frac_genus),
    names_to = "end",
    values_to = "frac"
  ) %>%
  group_by(depth, end) %>%
  mutate(
    # weighted mean
    w_mean = weighted.mean(frac, genus_count_assembly, na.rm = TRUE),
    
    # weighted variance → sd
    w_var  = sum(genus_count_assembly * (frac - w_mean)^2, na.rm = TRUE) /
      sum(genus_count_assembly, na.rm = TRUE),
    w_sd   = sqrt(w_var),
    
    # weighted Z-score
    z_w    = (frac - w_mean) / w_sd,
    
    # outlier definition
    outlier = abs(z_w) > 3
  ) %>%
  ungroup()

figMO <- plot_ly() %>%
  add_markers(
    data = plot_dat_MO,
    x = ~depth,
    y = ~frac,
    color = ~end,
    symbols = c("circle", "x"),
    symbol = ~outlier,
    size = ~ifelse(outlier, 12, 8),
    text = ~paste(
      "identifier:", identifier,
      "<br>depth:", depth,
      "<br>", end, ":", round(frac, 3),
      "<br>Weighted Z:", round(z_w, 2),
      "<br>Outlier:", outlier
    ),
    hoverinfo = "text"
  ) %>%
  layout(
    xaxis = list(title = "Depth"),
    yaxis = list(title = "Deamination fraction (genus)"),
    legend = list(title = list(text = "Series"))
  )

figMO




plot_dat_ST <- adnaST %>%
  select(depth, identifier, genus_count_assembly,
         deam3p_frac_genus, deam5p_frac_genus) %>%
  pivot_longer(
    cols = c(deam3p_frac_genus, deam5p_frac_genus),
    names_to = "end",
    values_to = "frac"
  ) %>%
  group_by(depth, end) %>%
  mutate(
    # weighted mean
    w_mean = weighted.mean(frac, genus_count_assembly, na.rm = TRUE),
    
    # weighted variance → sd
    w_var  = sum(genus_count_assembly * (frac - w_mean)^2, na.rm = TRUE) /
      sum(genus_count_assembly, na.rm = TRUE),
    w_sd   = sqrt(w_var),
    
    # weighted Z-score
    z_w    = (frac - w_mean) / w_sd,
    
    # outlier definition
    outlier = abs(z_w) > 3
  ) %>%
  ungroup()

figST <- plot_ly() %>%
  add_markers(
    data = plot_dat_ST,
    x = ~depth,
    y = ~frac,
    color = ~end,
    symbols = c("circle", "x"),
    symbol = ~outlier,
    size = ~ifelse(outlier, 12, 8),
    text = ~paste(
      "identifier:", identifier,
      "<br>depth:", depth,
      "<br>", end, ":", round(frac, 3),
      "<br>Weighted Z:", round(z_w, 2),
      "<br>Outlier:", outlier
    ),
    hoverinfo = "text"
  ) %>%
  layout(
    xaxis = list(title = "Depth"),
    yaxis = list(title = "Deamination fraction (genus)"),
    legend = list(title = list(text = "Series"))
  )

figST


# Write clean DF's
outlier_table_ST <- plot_dat_ST %>%
  filter(outlier) %>%
  arrange(depth, end, desc(z_w)) %>%
  select(
    depth,
    identifier,
    end,
    frac,
    genus_count_assembly,
    w_mean,
    w_sd,
    z_w
  )

outlier_table_MO <- plot_dat_MO %>%
  filter(outlier) %>%
  arrange(depth, end, desc(z_w)) %>%
  select(
    depth,
    identifier,
    end,
    frac,
    genus_count_assembly,
    w_mean,
    w_sd,
    z_w
  )

outlier_table_ST <- subset(outlier_table_ST, outlier_table_ST$z_w < 0)
outlier_table_MO <- subset(outlier_table_MO, outlier_table_MO$z_w < 0)


write.csv(outlier_table_MO, "C:/projects/repositories/zo_repo/output/outlier_table_MO.csv", row.names = FALSE)
write.csv(outlier_table_ST, "C:/projects/repositories/zo_repo/output/outlier_table_ST.csv", row.names = FALSE)

#####
dist_dat <- adnaST %>%
  select(depth, genus_count_assembly,
         deam3p_frac_genus, deam5p_frac_genus) %>%
  pivot_longer(
    cols = c(deam3p_frac_genus, deam5p_frac_genus),
    names_to = "end",
    values_to = "frac"
  )


ggplot(
  dist_dat,
  aes(
    x = frac,
    y = factor(depth),
    fill = end,
    weight = genus_count_assembly   # <- weights here
  )
) +
  geom_density_ridges(
    alpha = 0.7,
    scale = 3,
    rel_min_height = 0.01,
    color = "white"
  ) +
  labs(
    x = "Deamination fraction (genus)",
    y = "Depth",
    fill = "End",
    title = "Weighted distribution of deamination per depth"
  ) +
  theme_minimal()

library(ggridges)
dist_dat <- adnaMO %>%
  select(depth, genus_count_assembly,
         deam3p_frac_genus, deam5p_frac_genus) %>%
  pivot_longer(
    cols = c(deam3p_frac_genus, deam5p_frac_genus),
    names_to = "end",
    values_to = "frac"
  )

ggplot(
  dist_dat,
  aes(
    x = frac,
    y = factor(depth),
    fill = end,
    weight = genus_count_assembly   # <- weights here
  )
) +
  geom_density_ridges(
    alpha = 0.7,
    scale = 3,
    rel_min_height = 0.01,
    color = "white"
  ) +
  labs(
    x = "Deamination fraction (genus)",
    y = "Depth",
    fill = "End",
    title = "Weighted distribution of deamination per depth"
  ) +
  theme_minimal()

