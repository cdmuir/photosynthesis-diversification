source("r/header.R")

ephy <- read_rds("objects/ephy.rds")

tip_rate <- getTipRates(ephy, returnNetDiv = FALSE, statistic = "median")

# Check convergence and sample size
post <- as_draws(t(tip_rate$lambda))
diag_tbl <- summarize_draws(post, median, sd, rhat, ess_bulk, ess_tail)

diag_tbl |>
  arrange(desc(rhat))

hist(diag_tbl$rhat, breaks = 20)
