#!/usr/bin/env Rscript
# Read-only, independent verification; no fitting or source of analysis code.
options(warn=2)
self <- normalizePath(sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1]))
setwd(normalizePath(file.path(dirname(self),"../..")))
out <- "output/equal_validation_v1_20261016"
for (file in c("input_manifest.csv","output_manifest.csv")) {
  manifest <- read.csv(file.path(out,file))
  paths <- if(file=="input_manifest.csv") manifest$file else file.path(out,manifest$file)
  stopifnot(!anyDuplicated(paths),all(file.exists(paths)),
            all(file.info(paths)$size==manifest$bytes),
            identical(unname(vapply(paths,function(f)digest::digest(file=f,algo="sha256"),"")),
                      manifest$sha256))
  cat(file,":",nrow(manifest),"files verified\n")
}
z <- readRDS(file.path(out,"analysis.rds"))
p <- readRDS("simulation/results/arfc_select_v1/production_run/arfc_select_v1_run.rds")
b <- readRDS("simulation/results/arfc_legacy_benchmark_v1/production_run_v1/arfc_legacy_benchmark_v1_run_v1.rds")
libs <- list(fixed_select=c("fixed_fmcp","fixed_fscad"),
             asr_select=c("asr_mcp","asr_scad"),
             amr_select=c("amr_mcp","amr_scad"),
             arfc_select=c("asr_mcp","asr_scad","amr_mcp","amr_scad"))
cols <- c("task_id","method","validation_response_mse","response_prediction_mse")
raw <- rbind(p$candidate_metrics[cols],subset(b$metrics,method %in% libs$fixed_select)[cols])
stopifnot(nrow(z$choices)==1600L,nrow(z$metrics)==4000L,
          all(table(z$metrics$scenario,z$metrics$method)==100L))
for (i in seq_len(nrow(z$choices))) {
  a <- z$choices[i,]; methods <- libs[[a$selector]]
  d <- raw[raw$task_id==a$task_id & raw$method %in% methods,]
  minimum <- min(d$validation_response_mse)
  eligible <- d$method[abs(d$validation_response_mse-minimum)<=1e-12*(1+abs(minimum))]
  choice <- methods[methods %in% eligible][1]
  stopifnot(a$selected_component==choice,a$candidate_count==length(methods),!a$test_used)
  reported <- subset(z$metrics,task_id==a$task_id & method==a$selector)
  stopifnot(nrow(reported)==1L,
            reported$response_prediction_mse==d$response_prediction_mse[d$method==choice])
}
for (i in seq_len(nrow(z$summary))) {
  a <- z$summary[i,]
  values <- z$metrics[z$metrics$scenario==a$scenario & z$metrics$method==a$method,a$metric]
  values <- values[is.finite(values)]
  stopifnot(length(values)==a$n)
  if(length(values)) stopifnot(abs(mean(values)-a$mean)<1e-12)
}
for (i in seq_len(nrow(z$contrasts))) {
  a <- z$contrasts[i,]
  l <- z$metrics[z$metrics$scenario==a$scenario & z$metrics$method==a$left,]
  r <- z$metrics[z$metrics$scenario==a$scenario & z$metrics$method==a$right,]
  difference <- l[[a$metric]]-r[[a$metric]][match(l$task_id,r$task_id)]
  difference <- difference[is.finite(difference)]
  stopifnot(length(difference)==a$n)
  if(length(difference)) {
    estimate <- mean(difference); se <- sd(difference)/sqrt(length(difference))
    interval <- estimate+c(-1,1)*qt(.975,length(difference)-1)*se
    stopifnot(abs(estimate-a$mean)<1e-12,abs(se-a$se)<1e-12,
              max(abs(interval-c(a$lower,a$upper)))<1e-12)
  }
}
stopifnot(all(aggregate(count~scenario+selector,z$frequencies,sum)$count==100L))
cat("PASS: 1600 independent selections, 400 summary rows, 360 paired contrasts; all four cases retained.\n")
