## Reviewer 1, Point 1, option (i) (Botond, 2026-09-09): regenerate Figure 3 (birth-type re-parameterisation)
## WITH the three maternal BMI x ancestry interaction terms retained, so that Fig. 3 is literally the
## Fig. 2a fit re-parameterised. GLM column: complete-case model C (fit_birthtype(TRUE), n = 3,482).
## Bayesian column: the Fig. 2a joint model with imputation (imputalt.rData, n = 4,951, BMI x ancestry
## included), re-parameterised draw by draw: non-first vaginal = -b_first, first cesarean = b_CD + b_int,
## non-first cesarean = b_CD - b_first; every other coefficient is the same draw. This is exact for the
## linear predictor; only the (vague) priors are stated on the other parameterisation.
## Figure style = make_fig2.R panel B (two_model_forest, tm_double, 5600 x 4000). Also writes the full
## coefficient table and an old-vs-new comparison against the published Fig. 3 (no BMI x ancestry).
## Run from this folder: Rscript point1_fig3_withBMIxAnc.R > point1_fig3_withBMIxAnc.log 2>&1
Sys.setenv(REV2_ROOT = normalizePath("..")); source("../_common/setup.R")
pdf(NULL)
REV <- hapo_bayes_dir()

getB <- function(path) {                       # posterior draws, columns named by JointAI varname
  e <- new.env(); load(path, envir = e); m <- get(ls(e)[1], envir = e)
  cl <- m$coef_list$CPbinary; M <- do.call(rbind, lapply(m$MCMC, as.matrix))
  B <- M[, cl$coef, drop = FALSE]; colnames(B) <- cl$varname
  cat(basename(path), ": rows", nrow(m$data), " NA cells", sum(is.na(m$data)), " draws", nrow(B), " betas", ncol(B), "\n"); B
}
col <- function(B, p) { k <- grep(p, colnames(B)); stopifnot(length(k) == 1); B[, k] }
tpstr <- function(z) sprintf("%.2e", 2 * pnorm(-abs(mean(z) / sd(z))))   # as in make_fig2.R
ortxt <- function(e, l, u) paste0(round(exp(e), 2), " (", round(exp(l), 2), " - ", round(exp(u), 2), ")")

two_model_forest <- function(glm_model, blist, var_labels, xlim, ticks, file) {
  g <- broom::tidy(glm_model, conf.int = TRUE)
  stopifnot(length(blist) == nrow(g), length(var_labels) == nrow(g))
  best <- sapply(blist, mean); blow <- sapply(blist, quantile, .025); bupp <- sapply(blist, quantile, .975)
  disp <- data.frame(Variables = var_labels,
    `Complete case estimates` = ortxt(g$estimate, g$conf.low, g$conf.high),
    `Bayesian mean estimates` = ortxt(best, blow, bupp),
    `OR (CIs* 2.5%-97.5%)` = paste(rep(" ", 40), collapse = " "),
    `Complete case p` = sprintf("%.2e", g$p.value), `Imputed p` = sapply(blist, tpstr),
    check.names = FALSE, stringsAsFactors = FALSE)
  p <- forest(data = disp, est = list(exp(g$estimate), exp(best)), lower = list(exp(g$conf.low), exp(blow)),
              upper = list(exp(g$conf.high), exp(bupp)), ci_column = 4, ref_line = 1,
              arrow_lab = c("Lower risk", "Higher risk"), xlim = xlim, ticks_at = ticks, theme = tm_double)
  save_forest(p, file, width = 5600, height = 4000)
  data.frame(term = g$term, label = var_labels, glm_OR = exp(g$estimate), glm_lo = exp(g$conf.low), glm_hi = exp(g$conf.high),
             glm_p = g$p.value, bayes_OR = exp(best), bayes_lo = exp(blow), bayes_hi = exp(bupp),
             bayes_tailp = as.numeric(sapply(blist, tpstr)), stringsAsFactors = FALSE)
}
labels20 <- c("(Intercept)", "Maternal age at OGTT", "Maternal BMI at OGTT", "Ethnicity - Afro-Caribbean",
  "Ethnicity - Thai", "Ethnicity - Hispanic", "Maternal GDM", "Age of gestation at delivery", "Birthweight",
  "Neonatal sex (female)", "Maternal HbA1c at OGTT", "Neonatal head circumference", "Neonatal length",
  "HOMA2-IR at OGTT", "Maternal weight gain (till OGTT)", "Neonatal Cord PG concentration",
  "Non-first born, vaginal delivery", "First born, cesarean delivery", "Non-first born, cesarean delivery",
  "Cesarean delivery * Neonatal Cord PG concentration")
labels23 <- c(labels20, "Maternal BMI at OGTT * Ethnicity - Afro-Caribbean",
              "Maternal BMI at OGTT * Ethnicity - Thai", "Maternal BMI at OGTT * Ethnicity - Hispanic")

## ---- NEW Figure 3: model C (with BMI x ancestry) + Fig 2a joint model re-parameterised ----------------
mC <- fit_birthtype(TRUE); btd <- bt_data()$d
cat("\nModel C: n =", nobs(mC), " AUC =", auc_of(mC, btd), " deviance =", round(deviance(mC), 3), "\n")
mA <- fit_final(); cat("Fig 2a:  n =", nobs(mA), " AUC =", auc_of(mA, d_cc_s), " deviance =", round(deviance(mA), 3),
                       " max|fitted diff| =", signif(max(abs(fitted(mA) - fitted(mC))), 3), "\n")
BA <- getB(file.path(REV, "imputalt.rData"))
stopifnot(any(grepl("MaternalBMIatOGTT:Ethnicity", colnames(BA))))
bf <- col(BA, "^firsttimepregnantTRUE$"); bc <- col(BA, "^CesareanSection1$")
bi <- col(BA, "^firsttimepregnantTRUE:CesareanSection1$|^CesareanSection1:firsttimepregnantTRUE$")
blistC <- list(col(BA, "^\\(Intercept\\)$"), col(BA, "^MaternalageatOGTT$"), col(BA, "^MaternalBMIatOGTT$"),
  col(BA, "^EthnicityBlack$"), col(BA, "^EthnicityAsian$"), col(BA, "^EthnicityHispanic$"), col(BA, "^MaternalGDM"),
  col(BA, "^Ageofgestationatdelivery$"), col(BA, "^Birthweight$"), col(BA, "^Neonatalsexfemale2$"),
  col(BA, "^mHbA1cpercent$"), col(BA, "^bNHCMn$"), col(BA, "^bNLNGMn$"), col(BA, "^HOMA2IR$"),
  col(BA, "^WtMgain.till.OGTT$"), col(BA, "^bCordPGCmmol.L$"),
  -bf, bc + bi, bc - bf,
  col(BA, "^CesareanSection1:bCordPGCmmol.L$|^bCordPGCmmol.L:CesareanSection1$"),
  col(BA, "^MaternalBMIatOGTT:EthnicityBlack$"), col(BA, "^MaternalBMIatOGTT:EthnicityAsian$"), col(BA, "^MaternalBMIatOGTT:EthnicityHispanic$"))
stopifnot(length(blistC) == length(coef(mC)))
## sanity: GLM term order of model C must be the order assumed above
stopifnot(grepl("Non-first born, vaginal", names(coef(mC))[17]), grepl("First born, Cesarean", names(coef(mC))[18]),
          grepl("Non-first born, Cesarean", names(coef(mC))[19]), grepl("Cesarean Section", names(coef(mC))[20]),
          all(grepl("BMI", names(coef(mC))[21:23])))
newC <- two_model_forest(mC, blistC, labels23, xlim = c(0, 10), ticks = c(0.5, 1, 2, 4, 6, 8), file = "Figure_3_reparam_model_withBMIxAnc.jpg")
newC$glm_AUC <- auc_of(mC, btd)
write.csv(newC, "Figure_3_withBMIxAnc_coefficients.csv", row.names = FALSE)

## ---- published Figure 3 for comparison: model B (no BMI x ancestry) + b2b joint model --------------------
mB <- fit_birthtype(FALSE); cat("\nModel B (published Fig 3): n =", nobs(mB), " AUC =", auc_of(mB, btd), "\n")
BB <- getB(file.path(REV, "imputalt_osszevont_IMPUTED_b2b.rData"))
bfB <- col(BB, "^firsttimepregnantTRUE$"); bcB <- col(BB, "^CesareanSection1$")
biB <- col(BB, "^firsttimepregnantTRUE:CesareanSection1$|^CesareanSection1:firsttimepregnantTRUE$")
blistB <- list(col(BB, "^\\(Intercept\\)$"), col(BB, "^MaternalageatOGTT$"), col(BB, "^MaternalBMIatOGTT$"),
  col(BB, "^EthnicityBlack$"), col(BB, "^EthnicityAsian$"), col(BB, "^EthnicityHispanic$"), col(BB, "^MaternalGDM"),
  col(BB, "^Ageofgestationatdelivery$"), col(BB, "^Birthweight$"), col(BB, "^Neonatalsexfemale2$"),
  col(BB, "^mHbA1cpercent$"), col(BB, "^bNHCMn$"), col(BB, "^bNLNGMn$"), col(BB, "^HOMA2IR$"),
  col(BB, "^WtMgain.till.OGTT$"), col(BB, "^bCordPGCmmol.L$"),
  -bfB, bcB + biB, bcB - bfB,
  col(BB, "^CesareanSection1:bCordPGCmmol.L$|^bCordPGCmmol.L:CesareanSection1$"))
stopifnot(length(blistB) == length(coef(mB)))
oldB <- two_model_forest(mB, blistB, labels20, xlim = c(0, 10), ticks = c(0.5, 1, 2, 4, 6, 8), file = "_check_published_Fig3_regenerated.jpg")
oldB$glm_AUC <- auc_of(mB, btd)

## ---- old vs new, row by row ---------------------------------------------------------------------------------
f2 <- function(o, l, h) sprintf("%.2f (%.2f-%.2f)", o, l, h)
cmp <- merge(data.frame(label = oldB$label, old_GLM = f2(oldB$glm_OR, oldB$glm_lo, oldB$glm_hi), old_Bayes = f2(oldB$bayes_OR, oldB$bayes_lo, oldB$bayes_hi)),
             data.frame(label = newC$label, new_GLM = f2(newC$glm_OR, newC$glm_lo, newC$glm_hi), new_Bayes = f2(newC$bayes_OR, newC$bayes_lo, newC$bayes_hi)),
             by = "label", all = TRUE)
cmp <- cmp[match(labels23, cmp$label), ]
cmp$GLM_changed <- ifelse(is.na(cmp$old_GLM), "new row", ifelse(cmp$old_GLM == cmp$new_GLM, "", "CHANGED"))
cmp$Bayes_changed <- ifelse(is.na(cmp$old_Bayes), "new row", ifelse(cmp$old_Bayes == cmp$new_Bayes, "", "CHANGED"))
write.csv(cmp, "Figure_3_old_vs_new.csv", row.names = FALSE)
cat("\nOLD (published Fig 3, no BMI x ancestry) vs NEW (with BMI x ancestry), rounded as displayed:\n")
print(cmp, right = FALSE, row.names = FALSE)
cat(sprintf("\nAUC: published Fig 3 %.4f -> new Fig 3 %.4f (= Fig 2a %.4f)\n", auc_of(mB, btd), auc_of(mC, btd), auc_of(mA, d_cc_s)))
cat("DONE\n")
