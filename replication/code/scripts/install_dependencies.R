#!/usr/bin/env Rscript
# Explicit, opt-in installation into the archive-local library only.
stopifnot(file.exists("run.sh"),dir.exists("study"))
dir.create("library",showWarnings=FALSE)
lib <- normalizePath("library")
.libPaths(c(lib,.libPaths()))
options(repos=c(CRAN="https://cloud.r-project.org"))
if(!requireNamespace("remotes",quietly=TRUE))
  install.packages("remotes",lib=lib)
versions <- c(Rcpp="1.1.1-1.1",RcppArmadillo="15.2.7-1",fda="6.3.0",
              glmnet="5.0",flare="1.8",psych="2.6.3")
for(p in names(versions)) {
  same <- requireNamespace(p,quietly=TRUE) && packageVersion(p)==package_version(versions[[p]])
  if(!same) remotes::install_version(p,version=versions[[p]],lib=lib,
    upgrade="never",dependencies=NA,build_vignettes=FALSE)
}
reporting <- c("knitr","kableExtra","ggplot2","cowplot","png","rmarkdown","pkgdown","testthat","caret")
missing <- reporting[!vapply(reporting,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) install.packages(missing,lib=lib)
cat("Dependency setup finished. Install the three included source packages next.\n")
