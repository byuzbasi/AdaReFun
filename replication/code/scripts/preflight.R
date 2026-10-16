#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1L || !dir.exists(args[1])) stop("Supply the local library path.")
.libPaths(c(normalizePath(args[1]),.libPaths()))
expected <- c(FMCP="0.2.0",Rcpp="1.1.1.1.1",RcppArmadillo="15.2.7.1",
  PFLR="1.1.0",fda="6.3.0",glmnet="5.0",flare="1.8",psych="2.6.3")
observed <- vapply(names(expected),function(p) if(requireNamespace(p,quietly=TRUE))
  as.character(packageVersion(p)) else "MISSING",character(1))
tab <- data.frame(package=names(expected),expected=unname(expected),observed=unname(observed))
print(tab,row.names=FALSE)
if(any(expected!=observed)) stop("Frozen runtime mismatch. Restore recorded versions before production.")
cat("PASS: frozen numerical dependency versions.\n")
