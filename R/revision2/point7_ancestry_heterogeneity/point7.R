## Reviewer 1, round 2, Point 7: formal test of ancestry x cesarean delivery (and ancestry x first delivery)
## heterogeneity in the final complete-case model (n = 3,482), instead of four separate stratified fits.
## LRT and Wald tests of the interaction sets; ancestry-specific ORs from the interaction model (common covariates);
## Cochran's Q / I2 across the published stratified estimates (Table 2) as a secondary check; the stratified
## fits reproduced for verification. Run: Rscript point7.R > point7.log 2>&1.
Sys.setenv(REV2_ROOT = normalizePath("..")); source("../_common/setup.R")
m0 <- fit_final()
cat("\nFinal model: n =", nobs(m0), " AUC =", auc_of(m0, d_cc_s), "\n")
cat("Ethnicity levels:", paste(levels(d_cc_s$Ethnicity), collapse = " | "), "\n")
print(table(d_cc_s$Ethnicity))
INT_CD <- "Ethnicity`:`Cesarean Section"; INT_FB <- "Ethnicity`:`1st time pregnant"
m_cd <- fit_final(ints = c(int3, INT_CD)); m_fb <- fit_final(ints = c(int3, INT_FB)); m_both <- fit_final(ints = c(int3, INT_CD, INT_FB))
cat("\ncoefficient names, ancestry x CD model:\n"); print(names(coef(m_cd)))

lrt <- function(m1, m2, label) { a <- anova(m1, m2, test = "LRT")
  data.frame(test = label, df = a$Df[2], deviance_difference = a$Deviance[2], p_LRT = a$`Pr(>Chi)`[2],
             AIC_without = AIC(m1), AIC_with = AIC(m2)) }
tests <- rbind(lrt(m0, m_cd, "Ancestry x cesarean delivery added to the final model"),
               lrt(m0, m_fb, "Ancestry x first delivery added to the final model"),
               lrt(m0, m_both, "Both interaction sets added"),
               lrt(m_fb, m_both, "Ancestry x cesarean delivery, given ancestry x first delivery"),
               lrt(m_cd, m_both, "Ancestry x first delivery, given ancestry x cesarean delivery"))
cat("\nLikelihood-ratio tests:\n"); print(tests, digits = 4, row.names = FALSE)

wald <- function(m, idx) { b <- coef(m)[idx]; V <- vcov(m)[idx, idx, drop = FALSE]
  stat <- as.numeric(t(b) %*% solve(V) %*% b); c(chisq = stat, df = length(b), p = pchisq(stat, length(b), lower.tail = FALSE)) }
lev <- levels(d_cc_s$Ethnicity); ref <- lev[1]
int_name <- function(m, l, other) { nm <- names(coef(m)); nm[nm == paste0("Ethnicity", l, ":", other) | nm == paste0(other, ":Ethnicity", l)] }
idx_cd <- names(coef(m_cd)) %in% unlist(lapply(lev[-1], int_name, m = m_cd, other = "`Cesarean Section`1"))
idx_fb <- names(coef(m_fb)) %in% unlist(lapply(lev[-1], int_name, m = m_fb, other = "`1st time pregnant`TRUE"))
stopifnot(sum(idx_cd) == length(lev) - 1, sum(idx_fb) == length(lev) - 1)
wt <- rbind(`Ancestry x cesarean delivery (Wald)` = wald(m_cd, idx_cd), `Ancestry x first delivery (Wald)` = wald(m_fb, idx_fb))
cat("\nWald tests of the interaction coefficient sets:\n"); print(wt, digits = 4)
tests <- rbind(tests, data.frame(test = rownames(wt), df = wt[, "df"], deviance_difference = wt[, "chisq"], p_LRT = wt[, "p"], AIC_without = NA, AIC_with = NA))
write.csv(tests, "point7_heterogeneity_tests.csv", row.names = FALSE)

## ancestry-specific ORs from the interaction models (same covariates and interactions for every ancestry)
b_cd <- find_term(m_cd, "^`Cesarean Section`1$")
b_cdfirst <- find_term(m_cd, "Cesarean Section`1:`1st time pregnant`TRUE$|1st time pregnant`TRUE:`Cesarean Section`1$")
b_first <- find_term(m_fb, "^`1st time pregnant`TRUE$")
z_cd <- setNames(rep(0, length(coef(m_cd))), names(coef(m_cd))); z_fb <- setNames(rep(0, length(coef(m_fb))), names(coef(m_fb)))
anc <- do.call(rbind, lapply(lev, function(l) {
  L <- z_cd; L[b_cd] <- 1; if (l != ref) L[int_name(m_cd, l, "`Cesarean Section`1")] <- 1
  nf <- contrast(m_cd, L); L[b_cdfirst] <- 1; fb <- contrast(m_cd, L)
  M <- z_fb; M[b_first] <- 1; if (l != ref) M[int_name(m_fb, l, "`1st time pregnant`TRUE")] <- 1
  fd <- contrast(m_fb, M)
  data.frame(ancestry = l, n = sum(d_cc_s$Ethnicity == l),
             CD_vs_vaginal_within_nonfirstborns = fmt(nf[1:3]), p_CD_nonfirst = nf[4],
             CD_vs_vaginal_within_firstborns = fmt(fb[1:3]), p_CD_first = fb[4],
             first_vs_nonfirst_within_vaginal = fmt(fd[1:3]), p_first = fd[4]) }))
cat("\nAncestry-specific odds ratios from the interaction models (at mean cord glucose; common covariates):\n")
print(anc, digits = 3, row.names = FALSE); write.csv(anc, "point7_ancestry_specific_ORs.csv", row.names = FALSE)

## the four separate stratified fits (Table 2 / Supplementary Fig. 5-8), for verification
strat <- do.call(rbind, lapply(lev, function(l) {
  d <- d_cc %>% filter(Ethnicity == l); vars <- setdiff(magyarazo_bovitett_r2, "Ethnicity")
  ds <- standardize_predictors(d, vars); m <- fit_final(vars, int3[!grepl("Ethnicity", int3)], ds)
  o <- orci(m); cd <- o[o$term == "`Cesarean Section`1", ]; fi <- o[o$term == "`1st time pregnant`TRUE", ]
  cm <- birthtype_contrasts(m)   # delta-method contrasts within the stratum (CD within firstborns = main effect + CD x first)
  data.frame(ancestry = l, n = nobs(m), CD_OR = fmt(c(cd$OR, cd$lo, cd$hi)), p_CD = cd$p,
             CD_first_OR = fmt(cm[5, 1:3]), p_CD_first = cm[5, 4],
             first_OR = fmt(c(fi$OR, fi$lo, fi$hi)), p_first = fi$p, AUC = auc_of(m, ds)) }))
cat("\nSeparate stratified fits (check against Table 2):\n"); print(strat, digits = 3, row.names = FALSE)
write.csv(strat, "point7_stratified_check.csv", row.names = FALSE)

## Cochran's Q / I2 across the stratified estimates (log OR, SE from the 95% CI)
q_test <- function(o, label) { y <- log(o$or); se <- (log(o$hi) - log(o$lo)) / (2 * qnorm(0.975)); wgt <- 1 / se^2
  pooled <- sum(wgt * y) / sum(wgt); Q <- sum(wgt * (y - pooled)^2); df <- length(y) - 1
  data.frame(test = label, Q = Q, df = df, p = pchisq(Q, df, lower.tail = FALSE), I2_percent = max(0, (Q - df) / Q) * 100,
             pooled_OR = exp(pooled), pooled_lo = exp(pooled - qnorm(0.975) / sqrt(sum(wgt))), pooled_hi = exp(pooled + qnorm(0.975) / sqrt(sum(wgt)))) }
cq <- rbind(q_test(data.frame(or = c(1.92, 2.94, 1.45, 4.04), lo = c(0.97, 0.92, 0.67, 1.38), hi = c(3.75, 8.75, 2.98, 11.89)), "Cesarean delivery, Table 2 stratified ORs"),
            q_test(data.frame(or = c(0.30, 0.63, 0.82, 0.64), lo = c(0.16, 0.31, 0.33, 0.13), hi = c(0.56, 1.28, 1.91, 2.38)), "First delivery, Table 2 stratified ORs"))
cat("\nCochran's Q across the published stratified estimates:\n"); print(cq, digits = 4, row.names = FALSE)
write.csv(cq, "point7_cochranQ.csv", row.names = FALSE)
cat("\nDONE\n")
