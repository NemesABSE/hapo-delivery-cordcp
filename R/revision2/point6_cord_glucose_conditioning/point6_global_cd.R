## Reviewer 1, Point 6 follow-up (2026-09-17): "global" cesarean effect with cord glucose omitted, i.e. the final
## complete-case model WITHOUT the cesarean x first-delivery interaction (the parameterisation of the Results
## sentence "OR 2.91, 95% CI 2.13 to 3.97"), refitted with cord plasma glucose (main effect + CD x cord PG) omitted.
## Run from this folder: Rscript point6_global_cd.R > point6_global_cd.log 2>&1
Sys.setenv(REV2_ROOT = normalizePath("..")); source("../_common/setup.R")
cat("int3 =", paste(int3, collapse = " | "), "\n")
no_cdxfirst <- int3[!grepl("1st time pregnant", int3, fixed = TRUE)]
fit_report <- function(vars, ints, label) {
  m <- fit_final(vars, ints)
  cat("\n[", label, "] n =", nobs(m), " AUC =", auc_of(m, d_cc_s), "\n  interactions:", paste(ints, collapse = " | "), "\n")
  o <- orci(m); print(o[grepl("Cesarean|1st time|CordPG", o$term), ], digits = 3, row.names = FALSE)
  cd <- o[o$term == "`Cesarean Section`1", ]
  data.frame(model = label, n = nobs(m), AUC = auc_of(m, d_cc_s), CD_OR = fmt(c(cd$OR, cd$lo, cd$hi)), CD_p = signif(cd$p, 3),
             first_OR = { f <- o[o$term == "`1st time pregnant`TRUE", ]; fmt(c(f$OR, f$lo, f$hi)) })
}
vars_noPG <- setdiff(magyarazo_bovitett_r2, "b_CordPGC_mmol.L")
res <- rbind(
  fit_report(magyarazo_bovitett_r2, int3, "Full final model (3 interactions), CD main effect = within non-firstborns"),
  fit_report(magyarazo_bovitett_r2, no_cdxfirst, "Final model without CD x first delivery (global CD effect), cord PG kept"),
  fit_report(magyarazo_bovitett_r2, no_cdxfirst[!grepl("CordPG", no_cdxfirst)], "Without CD x first delivery and without CD x cord PG, cord PG main effect kept"),
  fit_report(vars_noPG, no_cdxfirst[!grepl("CordPG", no_cdxfirst)], "Without CD x first delivery, cord PG omitted (main effect and CD x cord PG): GLOBAL glucose-free CD effect"),
  fit_report(vars_noPG, character(0), "No interactions at all, cord PG omitted"))
write.csv(res, "point6_global_cd.csv", row.names = FALSE)
cat("\nSUMMARY:\n"); print(res, right = FALSE, row.names = FALSE); cat("\nDONE\n")
