## Supplementary Table 7: maternal GDM and neonatal hyperinsulinemia as No. (%) with the chi-square test (continuity correction),
## Thai participants with missing versus known birth order (the two rows were printed as mean (SD)).
Sys.setenv(REV2_ROOT = normalizePath("..")); source("../_common/setup.R")
thai <- adatok_master %>% filter(Ethnicity == "Asian"); thai$miss <- is.na(thai$`1st time pregnant`)
cat("Thai n =", nrow(thai), " missing =", sum(thai$miss), " known =", sum(!thai$miss), "\n")
out <- list()
for (v in c("Maternal GDM", "CP_binary")) {
  x <- as.numeric(as.character(thai[[v]])); cat("\n", v, ": values", paste(sort(unique(x)), collapse = ","), " NA", sum(is.na(x)), "\n")
  tb <- table(miss = thai$miss, x); print(tb); ch <- chisq.test(tb, correct = TRUE); tt <- t.test(x ~ thai$miss, var.equal = TRUE)
  f <- function(m) sprintf("%d (%.1f)", sum(x[thai$miss == m] == 1, na.rm = TRUE), 100 * mean(x[thai$miss == m] == 1, na.rm = TRUE))
  cat("missing:", f(TRUE), " known:", f(FALSE), " mean(SD) missing:", sprintf("%.2f (%.2f)", mean(x[thai$miss], na.rm = TRUE), sd(x[thai$miss], na.rm = TRUE)),
      " chi-square P =", signif(ch$p.value, 3), " t test P =", signif(tt$p.value, 3), "\n")
  out[[v]] <- data.frame(variable = v, missing = f(TRUE), known = f(FALSE), n_missing = sum(!is.na(x[thai$miss])), n_known = sum(!is.na(x[!thai$miss])), p_chisq = ch$p.value, p_t = tt$p.value)
}
write.csv(do.call(rbind, out), "table7_binary.csv", row.names = FALSE)
## how the other categorical rows of the printed table were tested (hypertension printed P = 0.089)
h <- grep("ypertens", names(thai), value = TRUE); cat("\nhypertension columns:", h, "\n")
