options(repos = c(CRAN = "https://cran.r-project.org"))

ins = c("BiocManager", "shiny", "DT", "ggplot2", "data.table", "bslib","stringr", "tibble", "dplyr")
bio = c("Biostrings", "Rfastp")

install_regular_if_missing = function(p) {
  if (p %in% rownames(installed.packages()) == FALSE) {
    install.packages(p)
  }
}

install_bio_if_missing = function(p) {
  if (p %in% rownames(installed.packages()) == FALSE) {
    BiocManager::install(p)
  }
}

invisible(sapply(ins, install_regular_if_missing))
invisible(sapply(bio, install_bio_if_missing))