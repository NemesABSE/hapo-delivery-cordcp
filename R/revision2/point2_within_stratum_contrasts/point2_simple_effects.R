## Reviewer 1, Point 2: 2 x 2 simple-effect contrasts of the final model (Fig 2a fit = Fig 3 re-parameterised).
## Rows: cesarean vs vaginal within each birth-order stratum; non-first vs first within each delivery mode.
## Complete case: delta-method contrasts of the Fig 2a GLM (n = 3,482).
## Bayesian: the same linear combinations applied draw by draw to the JointAI joint model with imputation
## (imputalt.rData, n = 4,951, BMI x ancestry included = the current Fig 3 Bayesian column); if the separately
## fitted birth-type model (option B) is present as point1_birthtype_B.rData it is reported as well.
## Run from this folder: Rscript point2_simple_effects.R > point2_simple_effects.log 2>&1
Sys.setenv(REV2_ROOT = normalizePath("..")); source("../_common/setup.R")

## log-OR linear combinations in terms of b_first (first delivery), b_cd (cesarean), b_int (cd x first)
##   first vaginal = 0; non-first vaginal = -b_first; first CS = b_cd + b_int; non-first CS = b_cd - b_first
L <- list(
  `CS vs vaginal, firstborns`              = c(first = 0,  cd = 1, int = 1),
  `CS vs vaginal, non-firstborns`          = c(first = 0,  cd = 1, int = 0),
  `Non-first vs first, vaginal delivery`   = c(first = -1, cd = 0, int = 0),
  `Non-first vs first, cesarean delivery`  = c(first = -1, cd = 0, int = -1),
  `Non-first vaginal vs first vaginal (Fig 3 row)`   = c(first = -1, cd = 0, int = 0),
  `First CS vs first vaginal (Fig 3 row)`            = c(first = 0,  cd = 1, int = 1),
  `Non-first CS vs first vaginal (Fig 3 row)`        = c(first = -1, cd = 1, int = 0))

cat("\n[GLM] Fig 2a complete-case final model\n")
mA <- fit_final(); cat("n =", nobs(mA), " deviance =", round(deviance(mA), 3), " AUC =", auc_of(mA, d_cc_s), "\n")
nm <- names(coef(mA)); z <- setNames(rep(0, length(nm)), nm)
b_first <- find_term(mA, "^`1st time pregnant`TRUE$"); b_cd <- find_term(mA, "^`Cesarean Section`1$")
b_int <- find_term(mA, "Cesarean Section`1:`1st time pregnant`TRUE$|1st time pregnant`TRUE:`Cesarean Section`1$")
stopifnot(length(b_first) == 1, length(b_cd) == 1, length(b_int) == 1)
glm_tab <- t(sapply(L, function(w) { v <- z; v[b_first] <- w["first"]; v[b_cd] <- w["cd"]; v[b_int] <- w["int"]; contrast(mA, v) }))
colnames(glm_tab) <- c("OR", "lo", "hi", "p")
for (i in seq_len(nrow(glm_tab))) cat(sprintf("  %-48s %s  p = %.2e\n", rownames(glm_tab)[i], fmt(glm_tab[i, 1:3]), glm_tab[i, 4]))
## consistency: the within-CS birth-order contrast equals the ratio of the two Fig 3 cesarean rows
cat("  check: 4.39-type ratio (non-first CS / first CS) =", round(glm_tab["Non-first CS vs first vaginal (Fig 3 row)", "OR"] / glm_tab["First CS vs first vaginal (Fig 3 row)", "OR"], 4),
    " vs direct contrast", round(glm_tab["Non-first vs first, cesarean delivery", "OR"], 4), "\n")

## Bayesian draws
bayes_tab <- function(rdata, label) {
  if (!file.exists(rdata)) { cat("\n", label, ": file not found, skipped:", rdata, "\n"); return(NULL) }
  e <- new.env(); load(rdata, envir = e); m <- get(ls(e)[1], envir = e)
  cl <- m$coef_list$CPbinary; M <- do.call(rbind, lapply(m$MCMC, as.matrix))
  B <- M[, cl$coef, drop = FALSE]; colnames(B) <- cl$varname
  cat("\n", label, "\n  file:", rdata, "\n  data rows", nrow(m$data), " posterior draws", nrow(B), " imputed cells", sum(is.na(m$data)), "\n")
  g <- function(p) { k <- grep(p, colnames(B)); stopifnot(length(k) == 1); B[, k] }
  bf <- g("^firsttimepregnantTRUE$"); bc <- g("^CesareanSection1$"); bi <- g("CesareanSection1:firsttimepregnantTRUE|firsttimepregnantTRUE:CesareanSection1")
  out <- t(sapply(L, function(w) { zz <- w["first"] * bf + w["cd"] * bc + w["int"] * bi
    c(OR = exp(mean(zz)), lo = exp(quantile(zz, .025)), hi = exp(quantile(zz, .975)), p_neg = mean(zz < 0)) }))
  colnames(out) <- c("OR", "lo", "hi", "P(logOR<0)")
  for (i in seq_len(nrow(out))) cat(sprintf("  %-48s %s  P(OR<1) = %.3f\n", rownames(out)[i], fmt(out[i, 1:3]), out[i, 4]))
  out
}
bA <- bayes_tab(file.path(hapo_bayes_dir(), "imputalt.rData"), "[Bayes A] Fig 2a joint model with imputation (BMI x ancestry included; current Fig 3 Bayesian column)")
bB <- bayes_tab("point1_birthtype_B.rData", "[Bayes B] separately fitted birth-type joint model (option B), if present")

tab <- data.frame(contrast = rownames(glm_tab),
                  GLM_complete_case = apply(glm_tab[, 1:3], 1, fmt), GLM_p = signif(glm_tab[, 4], 3),
                  Bayes_joint_model_imputed = if (is.null(bA)) NA else apply(bA[, 1:3], 1, fmt),
                  Bayes_P_OR_below_1 = if (is.null(bA)) NA else round(bA[, 4], 3), stringsAsFactors = FALSE)
if (!is.null(bB)) { tab$Bayes_optionB <- apply(bB[, 1:3], 1, fmt) }
write.csv(tab, "point2_simple_effects.csv", row.names = FALSE)
cat("\nSUMMARY (OR, 95% CI / credible interval):\n"); print(tab, right = FALSE, row.names = FALSE)
cat("\nDONE\n")
