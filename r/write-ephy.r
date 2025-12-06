source("r/header.R")

getEventData(
  phy = "bamm/phy.tre",
  eventdata = "bamm/event_data.txt",
  burnin = 0.25,
  type = "diversification"
) |>
  write_rds("objects/ephy.rds")
