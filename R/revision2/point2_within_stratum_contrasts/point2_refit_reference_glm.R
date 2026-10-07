## Reviewer 1, Point 2: complete-case final model (Fig 3 parameterisation, BMI x ancestry retained)
## REFITTED with the four-level Birth type variable relevelled to two further reference categories:
##   ref = Non-first born, vaginal delivery  -> direct estimate of CS vs vaginal within non-firstborns
##   ref = First born, Cesarean section      -> direct estimate of non-first vs first within cesarean
## plus the published reference (First born, vaginal delivery) as the check. Same 3,482 rows, same terms.
## Run from this folder: Rscript point2_refit_reference_glm.R > point2_refit_reference_glm.log 2>&1
Sys.setenv(REV2_ROOT = normalizePath("..")); source("../_common/setup.R")

bt <- bt_data()
refs <- c("First born, vaginal delivery", "Non-first born, vaginal delivery", "First born, Cesarean section")
rows <- list()
for (ref in refs) {
  b <- bt; b$d$`Birth type` <- relevel(b$d$`Birth type`, ref = ref)
  m <- fit_birthtype(TRUE, bt = b)
  cat(sprintf("\n[ref = %s]  n = %d  deviance = %.3f  AUC = %.4f\n", ref, nobs(m), deviance(m), auc_of(m, b$d)))
  o <- orci(m); o <- o[grepl("^`Birth type`", o$term), ]
  o$term <- sub("^`Birth type`", "", o$term)
  for (i in seq_len(nrow(o))) cat(sprintf("  %-36s vs ref: %s  p = %.2e\n", o$term[i], fmt(c(o$OR[i], o$lo[i], o$hi[i])), o$p[i]))
  rows[[ref]] <- data.frame(reference = ref, level = o$term, OR = o$OR, lo = o$lo, hi = o$hi, p = o$p,
                            deviance = deviance(m), stringsAsFactors = FALSE)
}
tab <- do.call(rbind, rows); rownames(tab) <- NULL
write.csv(tab, "point2_refit_reference_glm.csv", row.names = FALSE)
cat("\nDeviance identical across references:", isTRUE(all.equal(min(tab$deviance), max(tab$deviance))), "\n")
cat("\n2 x 2 simple effects, each read directly from the refit whose reference defines it:\n")
pick <- function(ref, lvl) { r <- tab[tab$reference == ref & tab$level == lvl, ]; fmt(c(r$OR, r$lo, r$hi)) }
cat(sprintf("  %-40s %s\n", "CS vs vaginal, firstborns",             pick("First born, vaginal delivery",     "First born, Cesarean section")))
cat(sprintf("  %-40s %s\n", "CS vs vaginal, non-firstborns",         pick("Non-first born, vaginal delivery", "Non-first born, Cesarean section")))
cat(sprintf("  %-40s %s\n", "Non-first vs first, vaginal delivery",  pick("First born, vaginal delivery",     "Non-first born, vaginal delivery")))
cat(sprintf("  %-40s %s\n", "Non-first vs first, cesarean delivery", pick("First born, Cesarean section",     "Non-first born, Cesarean section")))
cat("\nDONE\n")
