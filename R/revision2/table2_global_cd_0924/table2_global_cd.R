## Table 2 rework (2026-09-24): for each sensitivity scenario the GLOBAL cesarean delivery effect (model without the
## cesarean delivery x first delivery interaction) next to the within-birth-order effects of the interaction model
## (non-firstborn = CD main effect; firstborn = CD + CD x first, delta method). First delivery reported from the same models.
## Scenarios: final model, alternative glycemic covariates (Supp Fig 1), cord PG omitted, 95th-percentile outcome (final-model
## structure), four ancestry strata (stratified fits as in section 11 of the pipeline), ancestry x CD interaction model.
## Rscript table2_global_cd.R > table2_global_cd.log 2>&1
Sys.setenv(REV2_ROOT = normalizePath("..")); source("../_common/setup.R")

drop_first_cd <- function(ints) ints[!grepl("1st time pregnant", ints)]
b_cd  <- "`Cesarean Section`1"; b_fd <- "`1st time pregnant`TRUE"
row_of <- function(scen, m_int, m_glob, extra_L = NULL) {
  cm <- birthtype_contrasts(m_int); o <- orci(m_glob); oi <- orci(m_int)
  g <- o[o$term == b_cd, ]; gf <- o[o$term == b_fd, ]; fi <- oi[oi$term == b_fd, ]
  lr <- anova(m_glob, m_int, test = "LRT")
  data.frame(scenario = scen, n = nobs(m_int), n_glob = nobs(m_glob),
    CD_global = fmt(c(g$OR, g$lo, g$hi)), p_CD_global = g$p,
    CD_nonfirst = fmt(cm["CD vs vaginal within non-firstborns", 1:3]), p_CD_nonfirst = cm["CD vs vaginal within non-firstborns", 4],
    CD_first = fmt(cm["CD vs vaginal within firstborns", 1:3]), p_CD_first = cm["CD vs vaginal within firstborns", 4],
    FD_global = fmt(c(gf$OR, gf$lo, gf$hi)), p_FD_global = gf$p,
    FD_int_vaginal = fmt(c(fi$OR, fi$lo, fi$hi)), p_FD_int_vaginal = fi$p,
    p_LRT_firstxCD = lr$`Pr(>Chi)`[2], stringsAsFactors = FALSE)
}
res <- list()

## 1. final model
m1 <- fit_final(); m1g <- fit_final(ints = drop_first_cd(int3)); stopifnot(nobs(m1) == 3482)
res[[1]] <- row_of("Final complete-case model", m1, m1g)

## 2. alternative glycemic covariates (pipeline section 10a)
sens_vars <- c(setdiff(magyarazo_bovitett_r2, c("HOMA2_IR", "Maternal GDM")),
               "m_FCP_nmol_l", "m_FPG_CLC_mmol_l", "m_OneHrPG_CLC_mmol_l", "m_TwoHrPG_CLC_mmol_l")
uj <- adatok %>% select(all_of(id_variable), m_FCP_nmol_l, m_FPG_CLC_mmol_l, m_OneHrPG_CLC_mmol_l, m_TwoHrPG_CLC_mmol_l)
ds <- rename_vars(adatok_master %>% left_join(uj, by = id_variable), lab_map); sens_vars <- rename_vec(sens_vars, lab_map)
ds <- standardize_predictors(na.omit(ds[, c(sens_vars, target_variable)]), sens_vars)
m2 <- glm(build_formula(target_variable, c(sens_vars, final_interactions)), family = binomial(), data = ds)
m2g <- glm(build_formula(target_variable, c(sens_vars, drop_first_cd(final_interactions))), family = binomial(), data = ds)
res[[2]] <- row_of("Alternative glycemic covariates (Supp Fig 1)", m2, m2g)

## 3. cord PG omitted (main effect and CD x cord PG), BMI x ancestry kept as in point6
v3 <- setdiff(magyarazo_bovitett_r2, "b_CordPGC_mmol.L"); i3 <- int3[!grepl("CordPGC", int3)]
m3 <- fit_final(vars = v3, ints = i3); m3g <- fit_final(vars = v3, ints = drop_first_cd(i3))
res[[3]] <- row_of("Cord plasma glucose omitted", m3, m3g)

## 4. 95th-percentile outcome, final-model structure
d95 <- adatok_cp95 %>% select(all_of(c(id_variable, magyarazo_bovitett_r2, "CP_binary95"))) %>% na.omit()
d95 <- standardize_predictors(d95, magyarazo_bovitett_r2)
m4 <- glm(build_formula("CP_binary95", c(magyarazo_bovitett_r2, int3)), family = binomial(), data = d95)
m4g <- glm(build_formula("CP_binary95", c(magyarazo_bovitett_r2, drop_first_cd(int3))), family = binomial(), data = d95)
res[[4]] <- row_of("95th-percentile outcome (final-model structure)", m4, m4g)

## 5. ancestry strata (pipeline section 11: full model, 2 final interactions, ancestry dropped)
for (l in c("White", "Black", "Hispanic", "Asian")) {
  d <- na.omit(adatok[adatok$Ethnicity == l, c(magyarazo_bovitett_r2, target_variable)])
  vars <- setdiff(magyarazo_bovitett_r2, "Ethnicity"); d <- standardize_predictors(d, vars)
  m <- glm(build_formula(target_variable, c(vars, final_interactions)), family = binomial(), data = d)
  mg <- glm(build_formula(target_variable, c(vars, drop_first_cd(final_interactions))), family = binomial(), data = d)
  res[[length(res) + 1]] <- row_of(paste("Ancestry stratum:", l), m, mg)
}

## 6. ancestry x CD interaction model (point7 m_cd), ancestry-specific CD effects with and without CD x first
INT_CD <- "Ethnicity`:`Cesarean Section"
mc <- fit_final(ints = c(int3, INT_CD)); mcg <- fit_final(ints = c(drop_first_cd(int3), INT_CD))
lev <- levels(d_cc_s$Ethnicity); ref <- lev[1]
iname <- function(m, l) { nm <- names(coef(m)); nm[nm == paste0("Ethnicity", l, ":", b_cd) | nm == paste0(b_cd, ":Ethnicity", l)] }
bint <- find_term(mc, "Cesarean Section`1:`1st time pregnant`TRUE$|1st time pregnant`TRUE:`Cesarean Section`1$")
anc <- do.call(rbind, lapply(lev, function(l) {
  z <- setNames(rep(0, length(coef(mc))), names(coef(mc))); L <- z; L[b_cd] <- 1; if (l != ref) L[iname(mc, l)] <- 1
  nf <- contrast(mc, L); L[bint] <- 1; fb <- contrast(mc, L)
  zg <- setNames(rep(0, length(coef(mcg))), names(coef(mcg))); G <- zg; G[b_cd] <- 1; if (l != ref) G[iname(mcg, l)] <- 1
  gl <- contrast(mcg, G)
  data.frame(ancestry = l, CD_global = fmt(gl[1:3]), p_CD_global = gl[4], CD_nonfirst = fmt(nf[1:3]), p_CD_nonfirst = nf[4],
             CD_first = fmt(fb[1:3]), p_CD_first = fb[4]) }))

out <- do.call(rbind, res)
options(width = 250)
cat("\n=== Table 2 effects: global (no CD x first) vs within birth order (interaction model) ===\n"); print(out, digits = 3, row.names = FALSE)
cat("\n=== Ancestry x CD interaction model: ancestry-specific CD effects (Wald CI, delta method) ===\n"); print(anc, digits = 3, row.names = FALSE)
write.csv(out, "table2_global_cd.csv", row.names = FALSE); write.csv(anc, "table2_global_cd_ancestry_interaction.csv", row.names = FALSE)
cat("DONE\n")

## 7. ancestry x first delivery interaction model (point7 m_fb): ancestry-specific FD, with and without CD x first
INT_FB <- "Ethnicity`:`1st time pregnant"
mf <- fit_final(ints = c(int3, INT_FB)); mfg <- fit_final(ints = c(drop_first_cd(int3), INT_FB))
fname <- function(m, l) { nm <- names(coef(m)); nm[nm == paste0("Ethnicity", l, ":", b_fd) | nm == paste0(b_fd, ":Ethnicity", l)] }
ancf <- do.call(rbind, lapply(lev, function(l) {
  z <- setNames(rep(0, length(coef(mf))), names(coef(mf))); L <- z; L[b_fd] <- 1; if (l != ref) L[fname(mf, l)] <- 1; v <- contrast(mf, L)
  zg <- setNames(rep(0, length(coef(mfg))), names(coef(mfg))); G <- zg; G[b_fd] <- 1; if (l != ref) G[fname(mfg, l)] <- 1; g <- contrast(mfg, G)
  data.frame(ancestry = l, FD_global = fmt(g[1:3]), p_FD_global = g[4], FD_vaginal = fmt(v[1:3]), p_FD_vaginal = v[4]) }))
cat("\n=== Ancestry x first delivery interaction model: ancestry-specific FD ===\n"); print(ancf, digits = 3, row.names = FALSE)
write.csv(ancf, "table2_global_fd_ancestry_interaction.csv", row.names = FALSE)
cat("DONE7\n")
