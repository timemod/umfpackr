#!/usr/bin/Rscript

if(require(covr)){ 
   cv <- covr::package_coverage("pkg")
   print(cv)
   print(subset(tally_coverage(cv), value == 0),row.names=FALSE)
} else {
    stop("covr not installed")
}
