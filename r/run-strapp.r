source("r/header.R")

ephy <- read_rds("objects/ephy.rds")
load("eCO2volution_C3total_Hist.RData")

df_traits <- SummaryStatZ |>
  rownames_to_column("species") |>
  select(species, Vcmax, Jmax) |>
  mutate(
    genus = str_extract(species, "^[^_]*"),
    genus = case_when(
      genus == "Lycopersicon" ~ "Solanum",
      genus == "Leptobalanus" ~ "Licania",
      genus == "Lagerstroemiaeciosa" ~ "Lagerstroemia",
      genus == "unknown" ~ NA_character_,
      TRUE ~ genus
    )
  ) |>
  filter(!is.na(genus)) |>
  summarize(
    Jmax = mean(Jmax, na.rm = TRUE),
    Vcmax = mean(Vcmax, na.rm = TRUE),
    .by = genus
  )

Vcmax <- df_traits$Vcmax |>
  set_names(df_traits$genus) |>
  na.omit()

Jmax <- df_traits$Jmax |>
  set_names(df_traits$genus) |>
  na.omit()

focal_tips <- unique(c(names(Jmax), names(Vcmax)))
post_rate1 <- getTipRates(subtreeBAMM(ephy, tips = focal_tips), statistic = "median")
post_rate2 <- getTipRates(subtreeBAMM(ephy, tips = focal_tips),
                          statistic = "median",
                          returnNetDiv = TRUE)

df_strapp <- full_join(
  post_rate1 |>
    extract2("lambda.avg") |>
    enframe(name = "genus", value = "lambda"),
  post_rate1 |>
    extract2("mu.avg") |>
    enframe(name = "genus", value = "mu"),
  by = join_by(genus)
) |>
  full_join(
    post_rate2 |>
      extract2("netdiv.avg") |>
      enframe(name = "genus", value = "netdiv"),
    by = join_by(genus)
  ) |>
  full_join(enframe(Vcmax, name = "genus", value = "Vcmax"), by = join_by(genus)) |>
  full_join(enframe(Jmax, name = "genus", value = "Jmax"), by = join_by(genus))

strapp_Vcmax_speciation <- traitDependentBAMM(
  ephy = subtreeBAMM(ephy, tips = names(Vcmax)),
  traits = log(Vcmax),
  reps = 1e3,
  rate = "speciation",
  return.full = FALSE,
  method = "pearson",
  logrates = TRUE,
  two.tailed = TRUE,
  traitorder = NA,
  nthreads = 10
)

strapp_Vcmax_diversification <- traitDependentBAMM(
  ephy = subtreeBAMM(ephy, tips = names(Vcmax)),
  traits = log(Vcmax),
  reps = 1e3,
  rate = "net diversification",
  return.full = FALSE,
  method = "pearson",
  logrates = TRUE,
  two.tailed = TRUE,
  traitorder = NA,
  nthreads = 10
)

strapp_Jmax_speciation <- traitDependentBAMM(
  ephy = subtreeBAMM(ephy, tips = names(Jmax)),
  traits = log(Jmax),
  reps = 1e3,
  rate = "speciation",
  return.full = FALSE,
  method = "pearson",
  logrates = TRUE,
  two.tailed = TRUE,
  traitorder = NA,
  nthreads = 10
)

strapp_Jmax_diversification <- traitDependentBAMM(
  ephy = subtreeBAMM(ephy, tips = names(Jmax)),
  traits = log(Jmax),
  reps = 1e3,
  rate = "speciation",
  return.full = FALSE,
  method = "pearson",
  logrates = TRUE,
  two.tailed = TRUE,
  traitorder = NA,
  nthreads = 10
)

# Write
write_rds(df_strapp, "objects/df_strapp.rds")
write_rds(strapp_Vcmax_speciation,
          "objects/strapp_Vcmax_speciation.rds")
write_rds(strapp_Vcmax_diversification,
          "objects/strapp_Vcmax_diversification.rds")
write_rds(strapp_Jmax_speciation,
          "objects/strapp_Jmax_speciation.rds")
write_rds(strapp_Jmax_diversification,
          "objects/strapp_Jmax_diversification.rds")
