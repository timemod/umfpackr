#!/usr/bin/Rscript

repo <- "https://cloud.r-project.org"
if (!require("pak")) {
  install.packages("pak", repos = repo)
}
pak::local_install_dev_deps("pkg", lib = .libPaths(), upgrade = FALSE)
