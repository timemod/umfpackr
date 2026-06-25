source("install_deps.R")

# Make sure all packages are up-to-date.
update.packages(
  lib.loc = Sys.getenv("R_LIBS_USER"),
  ask = FALSE
)
