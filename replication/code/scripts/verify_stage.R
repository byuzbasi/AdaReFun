#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1L)
files <- list.files(args[1L],pattern="manifest.*csv$",full.names=TRUE)
stopifnot(length(files)==1L)
m <- read.csv(files)
p <- file.path(args[1L],m$file)
stopifnot(all(file.exists(p)),!anyDuplicated(m$file),
  all(file.info(p)$size==m$bytes),all(unname(tools::md5sum(p))==m$md5))
cat("PASS: completed-stage manifest.\n")
