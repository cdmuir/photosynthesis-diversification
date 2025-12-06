source("r/header.R")

# Genus level tree from Dimitrov et al. (2023)
phy0 <- read.tree("bamm/Molecular_trees/Molecular_tree_140_210Ma_constraint.tree")
write.tree(phy0, file = "bamm/phy.tre")

# Sampling fractions per genus from Dimitrov et al. (2023)
cat("0.75\n", file = "bamm/sample_probs.txt", append = FALSE)
read_excel("bamm/41467_2023_43396_MOESM6_ESM.xlsx", sheet = "BAMM_sampling_fractions_global") |>
  filter(Genus %in% phy0$tip.label) |>
  mutate(Clade = Genus) |>
  select(Genus, Clade, Sampling_fraction) |>
  write_tsv("sample_probs.txt", append = TRUE, col_names = FALSE)

setBAMMpriors(phy0, total.taxa = 359208, outfile = NULL)

system("cd bamm; bamm -c myControlFile.txt")
