n_aso  <- 23231.0
n_vpa  <- 339497.0
n_syn  <- 7437.0
n_union <- 836595.0

expected   <- n_aso * n_vpa / n_union
enrichment <- n_syn / expected

cat("Expected:", round(expected, 1), "\n")
cat("Enrichment:", round(enrichment, 3), "x\n")
