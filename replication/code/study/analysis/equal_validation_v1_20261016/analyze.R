#!/usr/bin/env Rscript
# Saved-output analysis only. Never source a runner or fit a model.
options(stringsAsFactors = FALSE, warn = 2)
argv <- commandArgs(trailingOnly = FALSE)
self <- normalizePath(sub("^--file=", "", grep("^--file=", argv, value = TRUE)[1]))
code_dir <- dirname(self)
root <- normalizePath(file.path(code_dir, "../.."))
setwd(root)
out <- "output/equal_validation_v1_20261016"
stage <- paste0(out, "_staging")
stopifnot(!file.exists(out), !file.exists(stage), requireNamespace("digest", quietly=TRUE),
          requireNamespace("knitr", quietly=TRUE))
helper <- "simulation/50_analyze_arfc_legacy_benchmark_v1.R"
# Reuse only pure function definitions, never the historical script's top level.
for (expr in parse(helper)) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) && is.symbol(expr[[2]]) &&
      as.character(expr[[2]]) %in% c("verify_manifest", "mean_ci")) eval(expr)
}
dirs <- c("simulation/results/arfc_select_v1/production_run",
          "simulation/results/arfc_legacy_benchmark_v1/production_run_v1",
          "simulation/results/arfc_legacy_benchmark_v1/production_analysis_v1")
input_files <- c(self, file.path(code_dir,"PROTOCOL.txt"), helper,
                 "config/generated/v18/metrics_corrected_v18.csv")
for (d in dirs) {
  mf <- list.files(d, pattern="manifest.*csv$", full.names=TRUE)
  stopifnot(length(mf)==1L)
  m <- verify_manifest(d, mf)
  input_files <- c(input_files, mf, file.path(d, m$file))
}
input_files <- unique(input_files)
sha <- function(f) digest::digest(file=f, algo="sha256")
before <- data.frame(file=sub(paste0(root,"/"),"",input_files,fixed=TRUE),
                     bytes=as.numeric(file.info(input_files)$size),
                     sha256=vapply(input_files,sha,character(1)))
rownames(before) <- NULL
p <- readRDS(file.path(dirs[1],"arfc_select_v1_run.rds"))
b <- readRDS(file.path(dirs[2],"arfc_legacy_benchmark_v1_run_v1.rds"))
parent_mf <- file.path(dirs[1],"arfc_select_v1_manifest.csv")
stopifnot(p$mode=="production", b$mode=="production",
          identical(b$parent_signature,p$scientific_signature),
          identical(b$parent_manifest_md5,unname(tools::md5sum(parent_mf))),
          identical(p$tasks,b$tasks), nrow(p$tasks)==400L,
          all(table(p$tasks$scenario)==100L), all(p$tasks$n_train==150L),
          all(p$tasks$n_validation==150L), all(p$tasks$n_test==2000L),
          all(p$mass_checks$passed), !any(p$selection_outcomes$selection_used_test),
          all(b$diagnostics$data_provenance_max_difference<1e-12))
pc <- p$candidate_metrics
bc <- b$metrics[b$metrics$method %in% c("fixed_fmcp","fixed_fscad"),]
stopifnot(nrow(pc)==1600L,nrow(bc)==800L,
          all(with(pc, success & converged & refit_success & fit_usable & bic_selected &
                     !objective_increased)),
          all(with(bc, applicable & fit_success & converged & !test_used_for_tuning &
                     !validation_used_for_tuning)))
tp <- p$path_tuning
tb <- b$tuning[b$tuning$method %in% c("fixed_fmcp","fixed_fscad"),]
stopifnot(all(with(tb[tb$selected,], success & converged & refit_success)))
for (d in list(tp,tb)) {
  selected_col <- if ("bic_selected" %in% names(d)) "bic_selected" else "selected"
  for (g in split(d,interaction(d$task_id,d$method,drop=TRUE))) {
    chosen <- which(g[[selected_col]])
    ok <- with(g,success & converged & refit_success & is.finite(bic))
    stopifnot(length(chosen)==1L,ok[chosen],
              chosen==which(ok & g$bic==min(g$bic[ok]))[1L])
  }
}
ids <- c("task_id","scenario","case_label","replication",
         "seed_train","seed_validation","seed_test","method")
metrics <- c("response_prediction_mse","signal_prediction_mse","ise_total",
             "ise_null","ise_nonnull","true_zero_rate","false_zero_rate",
             "false_active_proportion","null_iou","boundary_error")
base <- rbind(pc[c(ids,"validation_response_mse",metrics)],
              bc[c(ids,"validation_response_mse",metrics)])
key <- function(d) paste(d$task_id,d$method,sep="::")
stopifnot(!anyDuplicated(key(base)),all(is.finite(base$validation_response_mse)),
          all(is.finite(base$response_prediction_mse)))
corrected <- read.csv(input_files[4],check.names=FALSE)
idx <- match(key(base),key(corrected))
stopifnot(!anyNA(idx),!anyDuplicated(key(corrected)))
for (nm in setdiff(metrics,c("ise_null","ise_nonnull"))) {
  stopifnot(isTRUE(all.equal(base[[nm]],corrected[[nm]][idx],tolerance=1e-12)))
}
base[c("ise_null","ise_nonnull")] <- corrected[idx,c("ise_null","ise_nonnull")]
stopifnot(max(abs(base$ise_total-base$ise_null-base$ise_nonnull))<1e-12)
for (nm in ids[1:7]) stopifnot(identical(base[[nm]],corrected[[nm]][idx]))
libraries <- list(fixed_select=c("fixed_fmcp","fixed_fscad"),
                  asr_select=c("asr_mcp","asr_scad"),
                  amr_select=c("amr_mcp","amr_scad"),
                  arfc_select=c("asr_mcp","asr_scad","amr_mcp","amr_scad"))
choose <- function(validation, priority) {
  stopifnot(identical(names(validation),c("task_id","method","validation_response_mse")))
  do.call(rbind,lapply(unique(validation$task_id),function(id) {
    d <- validation[validation$task_id==id,]
    d <- d[match(priority,d$method),]
    stopifnot(!anyNA(d),nrow(d)==length(priority))
    minimum <- min(d$validation_response_mse)
    tie <- d$validation_response_mse<=minimum+1e-12*(1+abs(minimum))
    data.frame(task_id=id,selected_component=d$method[which(tie)[1]],
               selected_validation_mse=d$validation_response_mse[which(tie)[1]],
               minimum_validation_mse=minimum,tie_count=sum(tie),
               candidate_count=length(priority),validation_n=150L,test_used=FALSE)
  }))
}
# Synthetic selector checks: tie priority, genuine minimum, row-order invariance.
toy <- data.frame(task_id="toy",method=c("a","b"),validation_response_mse=c(1,1+1e-13))
stopifnot(choose(toy,c("b","a"))$selected_component=="b")
toy$validation_response_mse <- c(1,0.9)
stopifnot(choose(toy,c("a","b"))$selected_component=="b",
          identical(choose(toy,c("a","b")),choose(toy[2:1,],c("a","b"))))
audits <- list(); selected <- list()
for (name in names(libraries)) {
  a <- choose(base[c("task_id","method","validation_response_mse")],libraries[[name]])
  a$selector <- name
  a$scenario <- p$tasks$scenario[match(a$task_id,p$tasks$task_id)]
  z <- base[match(paste(a$task_id,a$selected_component,sep="::"),key(base)),]
  z$selected_component <- z$method; z$method <- name
  audits[[name]] <- a; selected[[name]] <- z
}
audit <- do.call(rbind,audits)
old_choice <- p$selection_outcomes$selected_method[
  match(audits$arfc_select$task_id,p$selection_outcomes$task_id)]
stopifnot(identical(audits$arfc_select$selected_component,old_choice))
base$selected_component <- base$method
all_metrics <- rbind(base,do.call(rbind,selected)); rownames(all_metrics) <- NULL
stopifnot(nrow(all_metrics)==4000L,!anyDuplicated(key(all_metrics)))
labels <- c(fixed_fmcp="fMCP",fixed_fscad="fSCAD",fixed_select="Fixed-Select",
            asr_mcp="ASR-MCP",asr_scad="ASR-SCAD",asr_select="ASR-Select",
            amr_mcp="AMR-MCP",amr_scad="AMR-SCAD",amr_select="AMR-Select",
            arfc_select="ARFC-Select")
scenarios <- p$scenarios$scenario
summaries <- list(); contrasts <- list(); k <- 0L; j <- 0L
pairs <- data.frame(left=c("asr_mcp","asr_scad","amr_mcp","amr_scad",
                          "asr_select","amr_select",rep("arfc_select",3)),
                    right=c("fixed_fmcp","fixed_fscad","asr_mcp","asr_scad",
                            "fixed_select","asr_select","asr_select","amr_select","fixed_select"))
for (s in scenarios) {
  for (m in names(labels)) for (metric in metrics) {
    v <- all_metrics[all_metrics$scenario==s & all_metrics$method==m,metric]
    stopifnot(length(v)==100L,sum(is.finite(v)) %in% c(0L,100L))
    k <- k+1L; summaries[[k]] <- cbind(data.frame(scenario=s,method=m,metric=metric),
                                     as.data.frame(as.list(mean_ci(v))))
  }
  for (i in seq_len(nrow(pairs))) for (metric in metrics) {
    l <- all_metrics[all_metrics$scenario==s & all_metrics$method==pairs$left[i],]
    r <- all_metrics[all_metrics$scenario==s & all_metrics$method==pairs$right[i],]
    r <- r[match(l$task_id,r$task_id),]
    stopifnot(identical(l$task_id,r$task_id))
    delta <- l[[metric]]-r[[metric]]; ci <- mean_ci(delta)
    if (all(is.finite(delta)) && sd(delta)>1e-14) {
      stopifnot(max(abs(unname(t.test(delta)$conf.int)-ci[c("lower","upper")]))<1e-10)
    }
    j <- j+1L
    contrasts[[j]] <- cbind(data.frame(scenario=s,left=pairs$left[i],right=pairs$right[i],metric=metric),
                            as.data.frame(as.list(ci)))
  }
}
summary <- do.call(rbind,summaries); contrast <- do.call(rbind,contrasts)
frequency <- do.call(rbind,lapply(scenarios,function(s) do.call(rbind,lapply(names(libraries),function(m) {
  d <- audit[audit$scenario==s & audit$selector==m,]
  count <- as.integer(table(factor(d$selected_component,levels=libraries[[m]])))
  data.frame(scenario=s,selector=m,selected_component=libraries[[m]],
             count=count,percent=100*count/nrow(d))
}))))
ties <- aggregate(audit$tie_count>1,audit[c("scenario","selector")],sum)
names(ties)[3] <- "tied_tasks"
get_mean <- function(s,m,q) summary$mean[summary$scenario==s & summary$method==m & summary$metric==q]
s2 <- "case_ii_sparse"
table2 <- data.frame(Method=unname(labels),PMSE=vapply(names(labels),function(m)get_mean(s2,m,"response_prediction_mse"),0.0),
                     True_zero_pct=100*vapply(names(labels),function(m)get_mean(s2,m,"true_zero_rate"),0.0),
                     False_zero_pct=100*vapply(names(labels),function(m)get_mean(s2,m,"false_zero_rate"),0.0))
table_contrasts <- data.frame(Contrast=paste(labels[pairs$left],"minus",labels[pairs$right]))
for (metric in c("response_prediction_mse","true_zero_rate")) {
  d <- contrast[contrast$scenario==s2 & contrast$metric==metric,]
  mult <- if(metric=="true_zero_rate")100 else 1
  table_contrasts[[metric]] <- sprintf("%+.4g [%+.4g, %+.4g]",mult*d$mean,mult*d$lower,mult*d$upper)
}
names(table_contrasts)[2:3] <- c("PMSE difference [95% CI]","True-zero difference, pp [95% CI]")
safety <- do.call(rbind,lapply(seq_along(scenarios),function(i) {
  s <- scenarios[i]; ms <- names(libraries)
  data.frame(Case=paste("Case",c("I","II","III","IV")[i]),Method=labels[ms],
    PMSE=vapply(ms,function(m)get_mean(s,m,"response_prediction_mse"),0.0),
    False_active_pct=100*vapply(ms,function(m)get_mean(s,m,"false_active_proportion"),0.0),
    False_zero_pct=100*vapply(ms,function(m)get_mean(s,m,"false_zero_rate"),0.0))
}))
stopifnot(identical(before$sha256,unname(vapply(input_files,sha,character(1)))))
dir.create(stage,recursive=TRUE)
emit <- function(d,name) write.csv(d,file.path(stage,paste0(name,".csv")),row.names=FALSE,na="NA")
emit(before,"input_manifest"); emit(all_metrics,"task_metrics"); emit(audit,"selection_audit")
emit(summary,"method_summary"); emit(contrast,"paired_contrasts"); emit(frequency,"selection_frequency")
emit(ties,"tie_counts"); emit(table2,"case_ii"); emit(safety,"safety")
saveRDS(list(protocol=readLines(file.path(code_dir,"PROTOCOL.txt")),metrics=all_metrics,
            summary=summary,contrasts=contrast,choices=audit,frequencies=frequency,ties=ties,
            labels=labels,source_signatures=c(p$scientific_signature,b$scientific_signature)),
        file.path(stage,"analysis.rds"))
writeLines(capture.output(sessionInfo()),file.path(stage,"sessionInfo.txt"))
file.copy(file.path(code_dir,"PROTOCOL.txt"),file.path(stage,"PROTOCOL.txt"))
for (nm in c("case_ii","contrasts","safety")) {
  tab <- switch(nm,case_ii=table2,contrasts=table_contrasts,safety=safety)
  digits <- switch(nm,case_ii=c(0,7,2,2),contrasts=5,safety=c(0,0,7,2,2))
  writeLines(knitr::kable(tab,format="latex",booktabs=TRUE,digits=digits,row.names=FALSE),
             file.path(stage,paste0(nm,".tex")))
}
html_table <- function(d,digits=5) as.character(knitr::kable(d,format="html",digits=digits,row.names=FALSE))
html <- c('<!doctype html><html lang="tr"><meta charset="utf-8"><title>AdaReFun: Eşit doğrulama bütçesi</title>',
 '<style>body{font:16px/1.5 system-ui;max-width:1100px;margin:40px auto;padding:0 20px;color:#202632}table{border-collapse:collapse;width:100%;margin:20px 0}th,td{padding:8px 10px;border-bottom:1px solid #ddd;text-align:right}th:first-child,td:first-child{text-align:left}th{background:#eef2f5}h1{font-size:28px}h2{font-size:21px;margin-top:32px}.note{color:#526071}</style>',
 '<h1>AdaReFun: Eşit doğrulama bütçeli karşılaştırma</h1>',
 '<p>Mevcut 400 görev; her senaryoda 100 tekrar. Yeni model, simülasyon veya bootstrap çalıştırılmadı.</p>',
 '<p>Fixed-Select, ASR-Select ve AMR-Select iki aday; ARFC-Select dört aday arasından aynı 150 doğrulama gözlemini kullanarak seçim yapar. Hepsi aynı 2000 bağımsız test gözleminde değerlendirilir.</p>',
 '<h2>Case II: yerel sıfır bölgesi</h2>',html_table(table2,c(0,7,2,2)),
 '<h2>Ayrı bileşen karşılaştırmaları</h2>',html_table(table_contrasts),
 '<p class="note">Fark = soldaki yöntem eksi sağdaki yöntem. PMSE farkında negatif, doğru-sıfır farkında pozitif değer daha iyidir. Aralıklar görev-eşleştirilmiş, noktasal %95 t aralıklarıdır; yeni doğrulayıcı veya çoklu karşılaştırma kontrollü kanıt değildir.</p>',
 '<h2>Dört senaryoda seçici güvenlik özeti</h2>',html_table(safety,c(0,0,7,2,2)),
 '<p class="note">NA: ilgili senaryoda uygulanamaz; sıfır olarak doldurulmadı. Senaryolar birleştirilmedi.</p>',
 '<h2>Case II seçim sayıları</h2>',html_table(frequency[frequency$scenario==s2,c("selector","selected_component","count")],0),
 '<h2>Yorum sınırı</h2><p>Bu sonradan planlanmış analiz, doğrulama örneklemi farkını giderir. Aday sayısını ve hesaplama maliyetini eşitlemez. Sabit yöntemlerden ASR’ye geçiş ağırlıklandırma ile kriter normalizasyonunu birlikte değiştirir; normalizasyonun tek başına etkisi ölçülmez. Sıfır bölgeleri tahmin amaçlı seçim sonrasında değerlendirilir; seçimde gerçek katsayı veya test sonucu kullanılmaz.</p>',
 '<p>Makale v18, paket ve GitHub dosyaları değiştirilmedi. Tam sayısal kayıtlar ve bütünlük manifesti aynı klasördedir.</p></html>')
writeLines(html,file.path(stage,"results.html"))
stopifnot(identical(before$sha256,unname(vapply(input_files,sha,character(1)))))
payload <- list.files(stage,full.names=TRUE)
manifest <- data.frame(file=basename(payload),bytes=as.numeric(file.info(payload)$size),
                       sha256=unname(vapply(payload,sha,character(1))))
emit(manifest,"output_manifest")
stopifnot(file.rename(stage,out))
cat("PASS: 400 paired tasks; 2400 base records; 1600 choices; original ARFC choices reproduced.\n")
cat("No fits or random draws. Inputs preserved. Output:",out,"\n")
print(table2,row.names=FALSE); print(table_contrasts,row.names=FALSE)
