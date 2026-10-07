## Point 2, internal curiosity (Botond, 2026-09-10): the 2 x 2 simple effects of the Fig 3 complete-case model
## refitted SEPARATELY in male and female neonates (sex term dropped, all other terms of the final model kept,
## Birth type relevelled so that every contrast is a directly estimated coefficient). Not for the manuscript.
## Run from this folder: Rscript point2_refit_reference_by_sex.R > point2_refit_reference_by_sex.log 2>&1
Sys.setenv(REV2_ROOT = normalizePath("..")); source("../_common/setup.R")

sexvar <- grep("sex", magyarazo_bovitett_r2, value = TRUE, ignore.case = TRUE); stopifnot(length(sexvar) == 1)
cat("sex variable:", sexvar, "\n"); print(table(adatok_master[[sexvar]], useNA = "ifany"))
fit_bt_sex <- function(sub, ref) {
  bt <- bt_data(sub); bt$vars <- setdiff(bt$vars, sexvar); bt$d[[sexvar]] <- NULL
  bt$d$`Birth type` <- relevel(bt$d$`Birth type`, ref = ref)
  list(m = fit_birthtype(TRUE, bt = bt), d = bt$d)
}
pick <- function(o, lvl) { r <- o[o$term == paste0("`Birth type`", lvl), ]; c(r$OR, r$lo, r$hi, r$p) }
rows <- list()
for (sx in sort(unique(na.omit(adatok_master[[sexvar]])))) {
  sub <- adatok_master[!is.na(adatok_master[[sexvar]]) & adatok_master[[sexvar]] == sx, ]
  out <- list()
  for (ref in c("First born, vaginal delivery", "Non-first born, vaginal delivery", "First born, Cesarean section")) {
    f <- fit_bt_sex(sub, ref); out[[ref]] <- list(o = orci(f$m), n = nobs(f$m), auc = auc_of(f$m, f$d), ev = sum(f$d[[target_variable]] == 1), bt = table(f$d$`Birth type`))
  }
  cat(sprintf("\n[%s = %s]  n = %d  events = %d  AUC = %.4f\n", sexvar, as.character(sx), out[[1]]$n, out[[1]]$ev, out[[1]]$auc)); print(out[[1]]$bt)
  eff <- rbind(
    `CS vs vaginal, firstborns`             = pick(out[["First born, vaginal delivery"]]$o,     "First born, Cesarean section"),
    `CS vs vaginal, non-firstborns`         = pick(out[["Non-first born, vaginal delivery"]]$o, "Non-first born, Cesarean section"),
    `Non-first vs first, vaginal delivery`  = pick(out[["First born, vaginal delivery"]]$o,     "Non-first born, vaginal delivery"),
    `Non-first vs first, cesarean delivery` = pick(out[["First born, Cesarean section"]]$o,     "Non-first born, Cesarean section"))
  for (i in seq_len(nrow(eff))) cat(sprintf("  %-40s %s  p = %.2e\n", rownames(eff)[i], fmt(eff[i, 1:3]), eff[i, 4]))
  rows[[as.character(sx)]] <- data.frame(sex = as.character(sx), n = out[[1]]$n, events = out[[1]]$ev, contrast = rownames(eff),
                                         OR = eff[, 1], lo = eff[, 2], hi = eff[, 3], p = eff[, 4], stringsAsFactors = FALSE)
}
tab <- do.call(rbind, rows); rownames(tab) <- NULL
write.csv(tab, "point2_refit_reference_by_sex.csv", row.names = FALSE)
cat("\nDONE\n")
