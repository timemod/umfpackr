#!/usr/bin/Rscript

options(repos = c(CRAN = "https://cloud.r-project.org"))

user_lib_dir <- Sys.getenv("R_LIBS_USER")
if (!dir.exists(user_lib_dir)) {
  dir.create(user_lib_dir, recursive = TRUE)
}

cat("library-pad:\n")
print(.libPaths())

# Install 'pak' if needed
if (!require("pak", character.only = TRUE, quietly = TRUE)) {
  install.packages("pak", dependencies = FALSE)
} else {
  unloadNamespace("pak")
}

pak::local_install_dev_deps("pkg", lib = .libPaths(), upgrade = FALSE)
