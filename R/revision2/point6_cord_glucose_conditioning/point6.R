## Reviewer 1, Point 6: cord plasma glucose is post-exposure. Refit the final model with cord glucose
## (main effect + CD x cord PG) omitted, and with gestational age / birthweight omitted.
## Run from this folder: Rscript point6.R > point6.log 2>&1
Sys.setenv(REV2_ROOT = normalizePath("..")); source("../_common/setup.R")

drop_fit <- function(drop, label) {
  vars <- setdiff(magyarazo_bovitett_r2, drop)
  ints <- int3; for (dv in drop) ints <- ints[!grepl(dv, ints, fixed = TRUE)]
  m <- fit_final(vars, ints)
  cat("\n[", label, "] n =", nobs(m), " AUC =", auc_of(m, d_cc_s), " terms dropped:", paste(drop, collapse = ", "), "\n")
  o <- orci(m); print(o[grepl("Cesarean|1st time", o$term), ], digits = 3, row.names = FALSE)
  cm <- birthtype_contrasts(m); print_contrasts(cm)
  data.frame(model = label, n = nobs(m), AUC = auc_of(m, d_cc_s),
             CD_main_effect_within_nonfirst = fmt(cm[4, 1:3]), CD_within_first = fmt(cm[5, 1:3]),
             first_CD_vs_first_vaginal = fmt(cm[2, 1:3]), nonfirst_CD_vs_first_vaginal = fmt(cm[3, 1:3]),
             nonfirst_vaginal_vs_first_vaginal = fmt(cm[1, 1:3]))
}
p6 <- rbind(
  drop_fit(character(0), "Full final model (reference)"),
  drop_fit("b_CordPGC_mmol.L", "Cord plasma glucose omitted (main effect and CD x cord PG)"),
  drop_fit("Age of gestation at delivery", "Gestational age omitted"),
  drop_fit("Birthweight", "Birthweight omitted"),
  drop_fit(c("b_CordPGC_mmol.L", "Age of gestation at delivery", "Birthweight"), "Cord PG, gestational age and birthweight all omitted"))
write.csv(p6, "point6_post_exposure_variables_omitted.csv", row.names = FALSE)
cat("\nSUMMARY (same 3,482 complete cases in every model):\n"); print(p6, right = FALSE, row.names = FALSE)

## unadjusted and glucose-only contrasts for context
cat("\nContext: cord glucose by mode of delivery in the complete-case set\n")
print(with(d_cc, tapply(b_CordPGC_mmol.L, `Cesarean Section`, function(x) sprintf("median %.2f (IQR %.2f to %.2f)", median(x), quantile(x, .25), quantile(x, .75)))))
cat("\nDONE\n")
