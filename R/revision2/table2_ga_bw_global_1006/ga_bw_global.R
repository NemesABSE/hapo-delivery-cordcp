## 2026-10-06: GLOBAL cesarean delivery and first delivery effects (model without the cesarean delivery x first delivery
## interaction) with gestational age or birthweight omitted, to complete Table 2 / Supp Table 5 in the same three views
## as the other scenarios (global | within non-firstborns | within firstborns). Same helper as table2_global_cd_0924.
## Rscript ga_bw_global.R > ga_bw_global.log 2>&1
Sys.setenv(REV2_ROOT = normalizePath("..")); source("../_common/setup.R")
drop_first_cd <- function(ints) ints[!grepl("1st time pregnant", ints)]
b_cd <- "`Cesarean Section`1"; b_fd <- "`1st time pregnant`TRUE"
row_of <- function(scen, v) {
  m <- fit_final(vars = v); mg <- fit_final(vars = v, ints = drop_first_cd(int3)); o <- orci(mg); cm <- birthtype_contrasts(m)
  g <- o[o$term == b_cd, ]; gf <- o[o$term == b_fd, ]
  data.frame(scenario = scen, n = nobs(mg), CD_global = fmt(c(g$OR, g$lo, g$hi)), p_CD_global = g$p,
    FD_global = fmt(c(gf$OR, gf$lo, gf$hi)), p_FD_global = gf$p,
    CD_nonfirst = fmt(cm["CD vs vaginal within non-firstborns", 1:3]), CD_first = fmt(cm["CD vs vaginal within firstborns", 1:3]),
    AUC_global = auc_of(mg, d_cc_s), stringsAsFactors = FALSE) }
out <- rbind(row_of("Final complete-case model", magyarazo_bovitett_r2),
             row_of("Gestational age omitted", setdiff(magyarazo_bovitett_r2, "Age of gestation at delivery")),
             row_of("Birthweight omitted", setdiff(magyarazo_bovitett_r2, "Birthweight")))
print(out, row.names = FALSE); write.csv(out, "ga_bw_global.csv", row.names = FALSE); cat("DONE\n")
