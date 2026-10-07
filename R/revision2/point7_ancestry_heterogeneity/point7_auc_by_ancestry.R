## Point 7 follow-up 2026-09-18: discrimination (AUC) within each ancestry, pooled models (final model and the
## ancestry-interaction models behind the blue squares of the forest) versus the four separate stratified fits (red).
## Same complete cases (n = 3,482); within-ancestry AUC = ROC of the fitted values restricted to that ancestry.
Sys.setenv(REV2_ROOT = normalizePath("..")); source("../_common/setup.R")
INT_CD <- "Ethnicity`:`Cesarean Section"; INT_FB <- "Ethnicity`:`1st time pregnant"
mods <- list(final = fit_final(), plus_anc_x_CD = fit_final(ints = c(int3, INT_CD)),
             plus_anc_x_first = fit_final(ints = c(int3, INT_FB)), plus_both = fit_final(ints = c(int3, INT_CD, INT_FB)))
auc <- function(y, p) as.numeric(pROC::auc(pROC::roc(y, p, quiet = TRUE)))
y <- d_cc_s[[target_variable]]; lev <- levels(d_cc_s$Ethnicity)
res <- do.call(rbind, lapply(c("ALL", lev), function(l) {
  idx <- if (l == "ALL") rep(TRUE, nrow(d_cc_s)) else d_cc_s$Ethnicity == l
  row <- data.frame(ancestry = l, n = sum(idx), events = sum(as.numeric(as.character(y[idx])) == 1))
  for (nm in names(mods)) row[[paste0("AUC_", nm)]] <- auc(y[idx], fitted(mods[[nm]])[idx])
  if (l != "ALL") { d <- d_cc %>% filter(Ethnicity == l); vars <- setdiff(magyarazo_bovitett_r2, "Ethnicity")
    ds <- standardize_predictors(d, vars); m <- fit_final(vars, int3[!grepl("Ethnicity", int3)], ds)
    row$AUC_separate_stratified_fit <- auc(ds[[target_variable]], fitted(m)); row$n_params_stratified <- length(coef(m))
  } else { row$AUC_separate_stratified_fit <- NA; row$n_params_stratified <- NA }
  row }))
cat("parameters: ", paste(names(mods), sapply(mods, function(m) length(coef(m))), collapse = " | "), "\n")
options(width = 220); print(res, digits = 4, row.names = FALSE); write.csv(res, "point7_auc_by_ancestry.csv", row.names = FALSE)
cat("\nDONE\n")
