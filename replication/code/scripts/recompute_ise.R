#!/usr/bin/env Rscript
# Saved-result correction only: no model fits, random draws, or seed changes.
args <- commandArgs(trailingOnly = TRUE)
destination <- if (length(args)) args[1L] else "config/generated/v18"
if (dir.exists(destination)) stop("Refusing to overwrite: ", destination)
source_dir <- "simulation/results/arfc_legacy_benchmark_v1/production_analysis_v1"
manifest_path <- file.path(source_dir, "arfc_legacy_benchmark_v1_analysis_manifest_v1.csv")
m <- read.csv(manifest_path)
paths <- file.path(source_dir, m$file)
stopifnot(all(file.exists(paths)), all(file.info(paths)$size == m$bytes),
          all(unname(tools::md5sum(paths)) == m$md5))
a <- readRDS(file.path(source_dir, "arfc_legacy_benchmark_v1_analysis_v1.rds"))
old <- a$combined_metrics
stopifnot(nrow(old) == 4800L, !anyDuplicated(old[c("task_id", "method")]),
          all(a$validation$passed))
new <- old
valid <- old$applicable & old$fit_success & is.finite(old$ise_total)
# The stored total and null errors use the same linear quadrature rule.
# Subtraction evaluates its complementary mask exactly, without reconstructing
# classical curves from their lower-resolution plotting grid.
new$ise_null[valid & is.na(old$ise_null)] <- 0
new$ise_nonnull[valid] <- old$ise_total[valid] - new$ise_null[valid]
stopifnot(all(new$ise_nonnull[valid] >= -1e-12),
          max(abs(new$ise_total[valid] - new$ise_null[valid] -
                  new$ise_nonnull[valid])) < 1e-12,
          identical(old[setdiff(names(old), c("ise_null", "ise_nonnull"))],
                    new[setdiff(names(new), c("ise_null", "ise_nonnull"))]))
# Independent coefficient-level check for all 1,600 adaptive candidate fits.
p <- readRDS("simulation/results/arfc_select_v1/production_run/arfc_select_v1_run.rds")
t <- seq(0, 1, length.out = p$config$support_grid_points)
g <- diff(t); w <- c(g[1]/2, (head(g,-1)+tail(g,-1))/2, tail(g,1)/2)
knots <- c(rep(0,4), seq(0,1,length.out=33)[2:32], rep(1,4))
B <- splines::splineDesign(knots, t, ord=4, outer.ok=TRUE)
truth <- function(s) switch(s,
  case_i_null = numeric(length(t)),
  case_ii_sparse = ifelse(t<=.3, 2*(1-t)*sin(2*pi*(t+.2)),
                         ifelse(t>=.7, 2*t*sin(2*pi*(t-.2)), 0)),
  case_iii_crossing = 7*t^3+2*sin(4*pi*t+.2),
  case_iv_away = 4*sqrt(t)+exp(t^2)*cos(3*pi*t))
cfs <- p$selected_coefficients
errors <- numeric(nrow(cfs))
for (i in seq_len(nrow(cfs))) {
  v <- (drop(B %*% as.numeric(cfs[i,grep("^b_",names(cfs))])) - truth(cfs$scenario[i]))^2
  j <- match(paste(cfs$task_id[i],cfs$method[i]), paste(new$task_id,new$method))
  mask <- if(cfs$scenario[i]=="case_i_null") rep(TRUE,length(t)) else
    if(cfs$scenario[i]=="case_ii_sparse") t>=.3 & t<=.7 else rep(FALSE,length(t))
  errors[i] <- max(abs(c(sum(w*v),sum(w[mask]*v[mask]),sum(w[!mask]*v[!mask])) -
                       as.numeric(new[j,c("ise_total","ise_null","ise_nonnull")])))
}
stopifnot(max(errors)<1e-10)
dir.create(destination,recursive=TRUE)
write.csv(new,file.path(destination,"metrics_corrected_v18.csv"),row.names=FALSE)
change <- data.frame(task_id=old$task_id,method=old$method,scenario=old$scenario,
  old_ise_null=old$ise_null,new_ise_null=new$ise_null,
  old_ise_nonnull=old$ise_nonnull,new_ise_nonnull=new$ise_nonnull,
  delta=new$ise_nonnull-old$ise_nonnull)
write.csv(change,file.path(destination,"ise_change_audit_v18.csv"),row.names=FALSE)
case <- new[new$scenario=="case_ii_sparse",]
order <- a$config$reporting_order
avg <- do.call(rbind,lapply(order,function(id) {
  d <- case[case$method==id,]
  data.frame(Method=a$config$method_labels[[id]],
    PMSE=mean(d$response_prediction_mse),ISE0=mean(d$ise_null),
    ISE1=mean(d$ise_nonnull),TZ=mean(d$true_zero_rate),FZ=mean(d$false_zero_rate))
}))
tab <- data.frame(Method=avg$Method,PMSE=sprintf("%.6f",avg$PMSE),
  ISE0=sprintf("%.4f",avg$ISE0),ISE1=sprintf("%.4f",avg$ISE1),
  TZ=ifelse(is.finite(avg$TZ),sprintf("%.1f\\%%",100*avg$TZ),"--"),
  FZ=ifelse(is.finite(avg$FZ),sprintf("%.1f\\%%",100*avg$FZ),"--"))
latex <- knitr::kable(tab,"latex",booktabs=TRUE,escape=FALSE,
  col.names=c("Method","PMSE","$\\mathrm{ISE}_0$","$\\mathrm{ISE}_1$","True zero","False zero"),
  caption=paste0("Case-II performance on the common test samples. PMSE is response prediction MSE; ",
    "the two ISE components partition the integration grid. Exact-zero entries are inapplicable to nonsparse estimators."),
  label="arfc-legacy-caseii",align=c("l",rep("r",5)))
writeLines(as.character(latex),file.path(destination,"table_case_ii_benchmark_v18.tex"))
write.csv(avg,file.path(destination,"case_ii_summary_v18.csv"),row.names=FALSE)
writeLines(c("No fits rerun. Only regional ISE bookkeeping changed.",
  "Null and complementary masks partition the original quadrature nodes; empty regions contribute zero.",
  "PMSE, total ISE, selection, structural metrics, coefficients and seeds are unchanged.",
  sprintf("All 1600 adaptive coefficient reconstructions verified; maximum error %.3g.",max(errors)),
  sprintf("Complete 4800-row benchmark includes 200 inapplicable truth-aware Oracle rows.")),
  file.path(destination,"validation.txt"))
cat("PASS: corrected 4800 metric rows; verified 1600 coefficient reconstructions.\n")
