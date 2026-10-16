#!/usr/bin/env Rscript
# Run from a prepared work directory. Never replaces an existing data file.
fetch <- function(url,dest,md5) {
  if(!file.exists(dest)) {
    dir.create(dirname(dest),recursive=TRUE,showWarnings=FALSE)
    tmp <- tempfile(tmpdir=dirname(dest))
    on.exit(unlink(tmp),add=TRUE)
    download.file(url,tmp,mode="wb",quiet=FALSE)
    stopifnot(unname(tools::md5sum(tmp))==md5)
    stopifnot(file.rename(tmp,dest))
  }
  stopifnot(unname(tools::md5sum(dest))==md5)
}
fetch("https://zenodo.org/api/records/13122321/files/1976_NSlibrary_withmetadata.csv/content",
  "real_data/data/raw/neospectra_v3/1976_NSlibrary_withmetadata.csv",
  "8cb55e45af554a69bf2b818c9ff31492")
truck <- "PFLR-master/data/truck.rda"
if(!file.exists(truck)) {
  td <- tempfile("pflr-source-"); dir.create(td)
  tarfile <- file.path(td,"PFLR_1.1.0.tar.gz")
  urls <- c("https://cran.r-project.org/src/contrib/PFLR_1.1.0.tar.gz",
            "https://cran.r-project.org/src/contrib/Archive/PFLR/PFLR_1.1.0.tar.gz")
  ok <- FALSE
  for(url in urls) {
    status <- tryCatch(download.file(url,tarfile,mode="wb"),error=function(e) 1L)
    if(identical(status,0L)) { ok <- TRUE; break }
  }
  if(!ok) stop("Could not retrieve PFLR 1.1.0 from CRAN.")
  utils::untar(tarfile,files="PFLR/data/truck.rda",exdir=td)
  original <- file.path(td,"PFLR/data/truck.rda")
  stopifnot(unname(tools::md5sum(original))=="597a06ef330f8ff33dc38a419676940d")
  dir.create(dirname(truck),recursive=TRUE,showWarnings=FALSE)
  stopifnot(file.copy(original,truck,overwrite=FALSE))
}
stopifnot(unname(tools::md5sum(truck))=="597a06ef330f8ff33dc38a419676940d")
cat("PASS: unchanged raw source hashes.\n")

cat("Next: prepare-data and baselines, or import the six existing analysis inputs.\n")
