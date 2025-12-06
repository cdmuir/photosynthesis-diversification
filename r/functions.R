# Modified From ChatGPT
nearest_taxon <- function(dmat, taxon) {
  
  # distances from focal to all others
  dvec <- dmat[taxon, ]
  
  # remove self-distance
  dvec <- dvec[names(dvec) != taxon]
  
  # nearest neighbor
  nearest <- names(which.min(dvec))
  dist <- min(dvec)
  
  list(nearest_taxon = nearest, distance = dist)
}
