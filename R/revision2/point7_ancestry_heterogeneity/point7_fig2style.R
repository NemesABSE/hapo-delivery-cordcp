## Point 7: Fig. 2a-style forest plot of the final complete-case model extended with the cesarean delivery x ancestry
## interaction terms (all terms shown). Blue = published final model (GLM), red = the same model + CD x ancestry.
## Run: Rscript point7_fig2style.R > point7_fig2style.log 2>&1
Sys.setenv(REV2_ROOT = normalizePath("..")); source("../_common/setup.R")
pdf(NULL)
ortxt <- function(e, l, u) ifelse(is.na(e), "", paste0(round(exp(e), 2), " (", round(exp(l), 2), " - ", round(exp(u), 2), ")"))
m0 <- fit_final(); g0 <- broom::tidy(m0, conf.int = TRUE)
m1 <- fit_final(ints = c(int3, "Ethnicity`:`Cesarean Section")); g1 <- broom::tidy(m1, conf.int = TRUE)
cat("final n =", nobs(m0), "AUC", auc_of(m0, d_cc_s), "| + CD x ancestry n =", nobs(m1), "AUC", auc_of(m1, d_cc_s), "\n")
cat("LRT: "); print(anova(m0, m1, test = "LRT"))
lab <- c("(Intercept)", "Maternal age at OGTT", "Maternal BMI at OGTT", "Ethnicity - Afro-Caribbean",
  "Ethnicity - Thai", "Ethnicity - Hispanic", "1st time delivery", "Cesarean Delivery", "Maternal GDM",
  "Age of gestation at delivery", "Birthweight", "Neonatal sex (female)", "Maternal HbA1c at OGTT",
  "Neonatal head circumference", "Neonatal length", "HOMA2-IR at OGTT", "Maternal weight gain (till OGTT)",
  "Neonatal Cord PG concentration", "1st time delivery * Cesarean Delivery",
  "Cesarean Delivery * Neonatal Cord PG concentration", "Maternal BMI at OGTT * Ethnicity - Afro-Caribbean",
  "Maternal BMI at OGTT * Ethnicity - Thai", "Maternal BMI at OGTT * Ethnicity - Hispanic",
  "Cesarean Delivery * Ethnicity - Afro-Caribbean", "Cesarean Delivery * Ethnicity - Thai", "Cesarean Delivery * Ethnicity - Hispanic")
chk <- c("Intercept", "age", "BMI", "Black", "Asian", "Hispanic", "1st time", "Cesarean", "GDM", "gestation", "Birthweight",
         "sex", "HbA1c", "HCMn", "LNGMn", "HOMA2", "Wtgain|WtMgain", "CordPGC", "1st time.*Cesarean|Cesarean.*1st time",
         "CordPGC.*Cesarean|Cesarean.*CordPGC", "BMI.*Black|Black.*BMI", "BMI.*Asian|Asian.*BMI", "BMI.*Hispanic|Hispanic.*BMI",
         "Black.*Cesarean|Cesarean.*Black", "Asian.*Cesarean|Cesarean.*Asian", "Hispanic.*Cesarean|Cesarean.*Hispanic")
stopifnot(nrow(g1) == length(lab), all(mapply(function(p, t) grepl(p, t), chk, g1$term)), all(g0$term == g1$term[1:nrow(g0)]))
pad <- function(x, n) c(x, rep(NA, n - length(x)))
e0 <- pad(g0$estimate, nrow(g1)); l0 <- pad(g0$conf.low, nrow(g1)); u0 <- pad(g0$conf.high, nrow(g1)); p0 <- pad(g0$p.value, nrow(g1))
disp <- data.frame(Variables = lab, `Final model` = ortxt(e0, l0, u0), `Final model + CD * ancestry` = ortxt(g1$estimate, g1$conf.low, g1$conf.high),
  `OR (95% CI)` = paste(rep(" ", 40), collapse = " "), `Final model p` = ifelse(is.na(p0), "", sprintf("%.2e", p0)),
  `+ CD * ancestry p` = sprintf("%.2e", g1$p.value), check.names = FALSE, stringsAsFactors = FALSE)
tm_p7 <- forest_theme(base_size = 8, ci_Theight = 0.3, refline_col = "#63666A", ci_col = c("#377eb8", "#e41a1c"),
  footnote_col = "#636363", footnote_fontface = "italic", legend_name = "Model",
  legend_value = c("Final model (GLM)", "+ CD x ancestry (GLM)"), vertline_lty = c("dashed", "dotted"),
  vertline_col = c("#d6604d", "#bababa"), row_spacing = unit(4, "mm"), core = list(padding = unit(c(4, 3), "mm")))
p <- forest(data = disp, est = list(exp(e0), exp(g1$estimate)), lower = list(exp(l0), exp(g1$conf.low)),
            upper = list(exp(u0), exp(g1$conf.high)), ci_column = 4, ref_line = 1,
            arrow_lab = c("Lower risk", "Higher risk"), xlim = c(0, 5), ticks_at = c(0.5, 1, 2, 3),
            footnote = "Complete-case GLM, n = 3,482. Likelihood-ratio test of the three CD x ancestry terms: chi2 = 2.46, 3 df, P = 0.48. Reference categories as in Fig. 2a.",
            theme = tm_p7)
save_forest(p, "point7_fig2style_CDxAncestry.jpg", width = 5600, height = 4400)
write.csv(data.frame(term = g1$term, label = lab, final_OR = exp(e0), final_lo = exp(l0), final_hi = exp(u0), final_p = p0,
                     int_OR = exp(g1$estimate), int_lo = exp(g1$conf.low), int_hi = exp(g1$conf.high), int_p = g1$p.value),
          "point7_fig2style_coefficients.csv", row.names = FALSE)
print(disp[c(4:8, 19, 24:26), c(1, 2, 3, 6)], row.names = FALSE); cat("DONE\n")
