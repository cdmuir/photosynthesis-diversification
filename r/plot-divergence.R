# Plot divergence time estimated from "full" tree versus pruned tree
source("r/header.R")
library(castor)
load("eCO2volution_C3total_Hist.RData")

# Taxa used for Schweiger & Schweiger analysis
taxa <- tree3$Scenario.3$tip.label

all(taxa %in% phy$tip.label)
taxa <- taxa[(taxa %in% phy$tip.label)]

# full distance matrices for all tips in each tree
d1_full <- cophenetic(tree3$Scenario.3)
d2_full <- cophenetic(phy) # takes a minute

df1 <- map_dfr(taxa, \(.x) {
  tibble(
    taxon = .x,
    nearest_taxon_pruned = nearest_taxon(d1_full, .x)$nearest_taxon,
    nearest_taxon_full = nearest_taxon(d2_full, .x)$nearest_taxon,
    dist_pruned = nearest_taxon(d1_full, .x)$distance,
    dist_full = nearest_taxon(d2_full, .x)$distance,
    n_node2 = ifelse(
      taxon %in% phy$tip.label & nearest_taxon_pruned %in% phy$tip.label,
      get_pairwise_distances(phy, taxon, nearest_taxon_pruned, as_edge_counts =
                               TRUE) - 1,
      NA
    )
  )
})

# average bias and r2
median(100 * df1$dist_pruned / df1$dist_full)
cor(log10(df1$dist_pruned), log10(df1$dist_full)) ^ 2

# average nodes removed
median(df1$n_node2 - 1, na.rm = TRUE)

fig1_panelb = ggplot(df1, aes(x = dist_full, y = dist_pruned)) +
  geom_point() +
  geom_abline(slope = 1,
              intercept = 0,
              linetype = "dashed") +
  labs(x = expression(atop("Nearest taxon distance [Myr]",Zanne~et~al.~(2014)~tree)), 
       y = expression(atop("Nearest taxon distance [Myr]",Scweiger~and~Schweiger~(2024)~tree))) +
  coord_equal() +
  scale_x_log10(labels = label_number()) +
  scale_y_log10(labels = label_number()) +
  theme_cowplot()

write_rds(fig1_panelb, "objects/fig1_panelb.rds")
