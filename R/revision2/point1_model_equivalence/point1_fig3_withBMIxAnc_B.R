## Reviewer 1, Point 1, option B figure: GLM column = complete-case model C (birth type + CD x cord PG +
## BMI x ancestry, n = 3,482); Bayesian column = the SEPARATELY fitted JointAI joint model with the same
## terms (point1_fit_jointai_birthtype_withBMIxAnc.R; imputalt_birthtype_withBMIxAnc_IMPUTED.rData, n = 4,951).
## Also compares this Bayesian column with option A (Fig. 2a joint model re-parameterised) row by row.
## Run from this folder after the fit: Rscript point1_fig3_withBMIxAnc_B.R > point1_fig3_withBMIxAnc_B.log 2>&1
Sys.setenv(REV2_ROOT = normalizePath("..")); source("../_common/setup.R")
pdf(NULL)
REV <- hapo_bayes_dir()
getB <- function(path) {
  e <- new.env(); load(path, envir = e); m <- get(ls(e)[1], envir = e)
  cl <- m$coef_list$CPbinary; M <- do.call(rbind, lapply(m$MCMC, as.matrix))
  B <- M[, cl$coef, drop = FALSE]; colnames(B) <- cl$varname
  cat(basename(path), ": rows", nrow(m$data), " NA cells", sum(is.na(m$data)), " draws", nrow(B), " betas", ncol(B), "\n"); B
}
col <- function(B, p) { k <- grep(p, colnames(B)); if (length(k) != 1) stop("pattern ", p, " matched ", length(k), ": ", paste(colnames(B), collapse = " | ")); B[, k] }
tpstr <- function(z) sprintf("%.2e", 2 * pnorm(-abs(mean(z) / sd(z))))
ortxt <- function(e, l, u) paste0(round(exp(e), 2), " (", round(exp(l), 2), " - ", round(exp(u), 2), ")")
two_model_forest <- function(glm_model, blist, var_labels, xlim, ticks, file) {
  g <- broom::tidy(glm_model, conf.int = TRUE)
  stopifnot(length(blist) == nrow(g), length(var_labels) == nrow(g))
  best <- sapply(blist, mean); blow <- sapply(blist, quantile, .025); bupp <- sapply(blist, quantile, .975)
  disp <- data.frame(Variables = var_labels, `Complete case estimates` = ortxt(g$estimate, g$conf.low, g$conf.high),
    `Bayesian estimates` = ortxt(best, blow, bupp), `OR (CIs* 2.5%-97.5%)` = paste(rep(" ", 40), collapse = " "),
    `Complete case p` = sprintf("%.2e", g$p.value), `Imputed p` = sapply(blist, tpstr), check.names = FALSE, stringsAsFactors = FALSE)
  p <- forest(data = disp, est = list(exp(g$estimate), exp(best)), lower = list(exp(g$conf.low), exp(blow)),
              upper = list(exp(g$conf.high), exp(bupp)), ci_column = 4, ref_line = 1,
              arrow_lab = c("Lower risk", "Higher risk"), xlim = xlim, ticks_at = ticks, theme = tm_double)
  save_forest(p, file, width = 5600, height = 4000)
  data.frame(term = g$term, label = var_labels, glm_OR = exp(g$estimate), glm_lo = exp(g$conf.low), glm_hi = exp(g$conf.high),
             glm_p = g$p.value, bayes_OR = exp(best), bayes_lo = exp(blow), bayes_hi = exp(bupp),
             bayes_tailp = as.numeric(sapply(blist, tpstr)), stringsAsFactors = FALSE)
}
labels23 <- c("(Intercept)", "Maternal age at OGTT", "Maternal BMI at OGTT", "Ethnicity - Afro-Caribbean",
  "Ethnicity - Thai", "Ethnicity - Hispanic", "Maternal GDM", "Age of gestation at delivery", "Birthweight",
  "Neonatal sex (female)", "Maternal HbA1c at OGTT", "Neonatal head circumference", "Neonatal length",
  "HOMA2-IR at OGTT", "Maternal weight gain (till OGTT)", "Neonatal Cord PG concentration",
  "Non-first born, vaginal delivery", "First born, cesarean delivery", "Non-first born, cesarean delivery",
  "Cesarean delivery * Neonatal Cord PG concentration", "Maternal BMI at OGTT * Ethnicity - Afro-Caribbean",
  "Maternal BMI at OGTT * Ethnicity - Thai", "Maternal BMI at OGTT * Ethnicity - Hispanic")

mC <- fit_birthtype(TRUE); btd <- bt_data()$d
cat("Model C: n =", nobs(mC), " AUC =", auc_of(mC, btd), "\n")
BJ <- getB(file.path(REV, "imputalt_birthtype_withBMIxAnc_IMPUTED.rData"))
cat("JointAI varnames:\n"); print(colnames(BJ))
blist <- list(col(BJ, "^\\(Intercept\\)$"), col(BJ, "^MaternalageatOGTT$"), col(BJ, "^MaternalBMIatOGTT$"),
  col(BJ, "^EthnicityBlack$"), col(BJ, "^EthnicityAsian$"), col(BJ, "^EthnicityHispanic$"), col(BJ, "^MaternalGDM"),
  col(BJ, "^Ageofgestationatdelivery$"), col(BJ, "^Birthweight$"), col(BJ, "^Neonatalsexfemale2$"),
  col(BJ, "^mHbA1cpercent$"), col(BJ, "^bNHCMn$"), col(BJ, "^bNLNGMn$"), col(BJ, "^HOMA2IR$"),
  col(BJ, "^WtMgain.till.OGTT$"), col(BJ, "^bCordPGCmmol.L$"),
  col(BJ, "^BirthtypeNon-firstborn,vaginaldelivery$"), col(BJ, "^BirthtypeFirstborn,Cesareansection$"),
  col(BJ, "^BirthtypeNon-firstborn,Cesareansection$"),
  col(BJ, "CesareanSection1:bCordPGCmmol.L$|^bCordPGCmmol.L:CesareanSection1$"),
  col(BJ, "^MaternalBMIatOGTT:EthnicityBlack$|^EthnicityBlack:MaternalBMIatOGTT$"),
  col(BJ, "^MaternalBMIatOGTT:EthnicityAsian$|^EthnicityAsian:MaternalBMIatOGTT$"),
  col(BJ, "^MaternalBMIatOGTT:EthnicityHispanic$|^EthnicityHispanic:MaternalBMIatOGTT$"))
resB <- two_model_forest(mC, blist, labels23, xlim = c(0, 10), ticks = c(0.5, 1, 2, 4, 6, 8), file = "Figure_3_reparam_model_withBMIxAnc_jointB.jpg")
resB$glm_AUC <- auc_of(mC, btd)
write.csv(resB, "Figure_3_withBMIxAnc_jointB_coefficients.csv", row.names = FALSE)

## compare Bayesian columns: option A (Fig 2a joint model re-parameterised) vs option B (separate fit)
A <- read.csv("Figure_3_withBMIxAnc_coefficients.csv", stringsAsFactors = FALSE)
f2 <- function(o, l, h) sprintf("%.2f (%.2f-%.2f)", o, l, h)
cmp <- data.frame(label = labels23, GLM = f2(resB$glm_OR, resB$glm_lo, resB$glm_hi),
                  Bayes_A_reparam = f2(A$bayes_OR, A$bayes_lo, A$bayes_hi), Bayes_B_separate_fit = f2(resB$bayes_OR, resB$bayes_lo, resB$bayes_hi))
cmp$differs <- ifelse(cmp$Bayes_A_reparam == cmp$Bayes_B_separate_fit, "", "differs")
write.csv(cmp, "Figure_3_bayes_optionA_vs_optionB.csv", row.names = FALSE)
cat("\nBayesian column, option A (re-parameterised Fig 2a joint model) vs option B (separate joint fit):\n")
print(cmp, right = FALSE, row.names = FALSE)
cat("DONE\n")
