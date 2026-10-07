## Reviewer 1, Point 1: are Fig 2a and Fig 3 the same model? Refit Fig 3 with BMI x ancestry retained;
## report the Bayesian multilevel odds ratios numerically.
## Run from this folder: Rscript point1.R > point1.log 2>&1
Sys.setenv(REV2_ROOT = normalizePath("..")); source("../_common/setup.R")

cat("\n[A] Fig 2a model: separate terms + CD*cordPG + CD*first + Ethnicity*BMI (complete case)\n")
mA <- fit_final(); cat("n =", nobs(mA), " AUC =", auc_of(mA, d_cc_s), "\n")
oa <- orci(mA); print(oa[grepl("Cesarean|1st time|BMI", oa$term), ], digits = 3, row.names = FALSE)
cat("\n[A] multilevel contrasts DERIVED from the Fig 2a coefficients (delta method):\n"); cA <- birthtype_contrasts(mA); print_contrasts(cA)

cat("\n[B] Fig 3 model as published: Birth type variable, BMI x ancestry NOT retained\n")
mB <- fit_birthtype(FALSE); cat("n =", nobs(mB), " AUC =", auc_of(mB, bt_data()$d), "\n")
ob <- orci(mB); print(ob[grepl("Birth type|Cesarean", ob$term), ], digits = 3, row.names = FALSE)

cat("\n[C] Fig 3 parameterization WITH BMI x ancestry retained (algebraically identical to A)\n")
mC <- fit_birthtype(TRUE); cat("n =", nobs(mC), " AUC =", auc_of(mC, bt_data()$d), "\n")
oc <- orci(mC); print(oc[grepl("Birth type|Cesarean|BMI", oc$term), ], digits = 3, row.names = FALSE)

cat("\nEquivalence: deviance A", round(deviance(mA), 4), "| C", round(deviance(mC), 4), "| B", round(deviance(mB), 4),
    "; residual df A", mA$df.residual, "C", mC$df.residual, "B", mB$df.residual, "\n")
cat("max |fitted(A) - fitted(C)| =", signif(max(abs(fitted(mA) - fitted(mC))), 3), "\n")
cat("LRT for the 3 BMI x ancestry terms (B vs C): chi2 =", round(deviance(mB) - deviance(mC), 3),
    " df =", mB$df.residual - mC$df.residual, " p =", signif(pchisq(deviance(mB) - deviance(mC), mB$df.residual - mC$df.residual, lower.tail = FALSE), 3), "\n")

## Bayesian multilevel ORs, n = 4,951 (posterior mean, 2.5th to 97.5th percentiles)
bayes_contrasts <- function(rdata, label) {
  e <- new.env(); load(rdata, envir = e); m <- get(ls(e)[1], envir = e)
  cl <- m$coef_list$CPbinary; M <- do.call(rbind, lapply(m$MCMC, as.matrix))
  B <- M[, cl$coef, drop = FALSE]; colnames(B) <- cl$varname
  cat("\n", label, "\n  file:", rdata, "\n  data rows", nrow(m$data), " posterior draws", nrow(B), " imputed cells", sum(is.na(m$data)),
      " BMI x ancestry terms present:", any(grepl("MaternalBMIatOGTT:Ethnicity", colnames(B))), "\n")
  orq <- function(z) c(exp(mean(z)), exp(quantile(z, .025)), exp(quantile(z, .975)))
  g <- function(p) { k <- grep(p, colnames(B)); stopifnot(length(k) == 1); B[, k] }
  bf <- g("^firsttimepregnantTRUE$"); bc <- g("^CesareanSection1$"); bi <- g("CesareanSection1:firsttimepregnantTRUE|firsttimepregnantTRUE:CesareanSection1")
  out <- rbind(`Non-first born, vaginal (ref first vaginal)` = orq(-bf), `First born, cesarean (ref first vaginal)` = orq(bc + bi),
               `Non-first born, cesarean (ref first vaginal)` = orq(bc - bf), `CD vs vaginal within non-firstborns` = orq(bc),
               `CD vs vaginal within firstborns` = orq(bc + bi))
  for (i in seq_len(nrow(out))) cat(sprintf("  %-48s %s\n", rownames(out)[i], fmt(out[i, ])))
  out
}
bx <- bayes_contrasts(file.path(hapo_bayes_dir(), "imputalt.rData"), "[Bayes, WITH BMI x ancestry] Fig 2a joint model, seed 2020, 3 chains x 1060")
by <- bayes_contrasts(file.path(hapo_bayes_dir(), "imputalt_osszevont_IMPUTED_b2b.rData"), "[Bayes, WITHOUT BMI x ancestry] b2b joint model (published Fig 3 red column)")

tab <- data.frame(contrast = rownames(cA),
  GLM_Fig2a_derived = apply(cA[, 1:3], 1, fmt), GLM_Fig3_published_noBMIxAnc = NA, GLM_Fig3_with_BMIxAnc = NA,
  Bayes_with_BMIxAnc = apply(bx, 1, fmt), Bayes_without_BMIxAnc = apply(by, 1, fmt))
getb <- function(o, p) { r <- o[grepl(p, o$term), ]; fmt(c(r$OR, r$lo, r$hi)) }
for (o in list(list(ob, "GLM_Fig3_published_noBMIxAnc"), list(oc, "GLM_Fig3_with_BMIxAnc"))) {
  tab[1, o[[2]]] <- getb(o[[1]], "Non-first born, vaginal"); tab[2, o[[2]]] <- getb(o[[1]], "First born, Cesarean"); tab[3, o[[2]]] <- getb(o[[1]], "Non-first born, Cesarean")
}
write.csv(tab, "point1_multilevel_estimates.csv", row.names = FALSE)
cat("\nSUMMARY TABLE (OR, 95% CI or credible interval):\n"); print(tab, right = FALSE, row.names = FALSE)
cat("\nDONE\n")
