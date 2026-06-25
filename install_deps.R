#!/usr/bin/Rscript

options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!require("pak")) {
  install.packages("pak")
}
pak::local_install_dev_deps("pkg", lib = .libPaths(), upgrade = FALSE)
