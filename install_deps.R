#!/usr/bin/Rscript
repo <- "https://cloud.r-project.org"
if (!require(devtools)) {
  install.packages("devtools", repos = repo)
}
devtools::install_deps("pkg", dependencies = TRUE,
                       upgrade = "never", repos = repo)

