## Reviewer 1, Point 2 follow-up (Botond, 2026-09-11): the Fig. 2a complete-case final model refitted WITHOUT the
## first delivery x cesarean delivery interaction (the other two interactions, CD x cord PG and BMI x ancestry, kept),
## so that the cesarean main effect is a single birth-order-averaged contrast. Same 3,482 rows as the final model.
## Run from this folder: Rscript point2_fig2a_no_firstxCD_interaction.R > point2_fig2a_no_firstxCD_interaction.log 2>&1
Sys.setenv(REV2_ROOT = normalizePath("..")); source("../_common/setup.R")
cat("final-model interactions:", paste(int3, collapse = " | "), "\n")
int_no <- int3[!grepl("1st time pregnant", int3)]
stopifnot(length(int_no) == length(int3) - 1)
m_full <- fit_final(); m_no <- fit_final(ints = int_no)
stopifnot(nobs(m_full) == 3482, nobs(m_no) == 3482)
tab <- function(m) { g <- broom::tidy(m, conf.int = TRUE); g$OR <- exp(g$estimate); g$lo <- exp(g$conf.low); g$hi <- exp(g$conf.high); g }
tn <- tab(m_no); tf <- tab(m_full)
cat(sprintf("\nFull final model:        n = %d  deviance = %.3f  AIC = %.2f  AUC = %.4f\n", nobs(m_full), deviance(m_full), AIC(m_full), auc_of(m_full, d_cc_s)))
cat(sprintf("Without first x CD:      n = %d  deviance = %.3f  AIC = %.2f  AUC = %.4f\n", nobs(m_no), deviance(m_no), AIC(m_no), auc_of(m_no, d_cc_s)))
lrt <- anova(m_no, m_full, test = "LRT"); cat(sprintf("LRT for the interaction: chi-square %.2f on %d df, P = %.4f\n", lrt$Deviance[2], lrt$Df[2], lrt$`Pr(>Chi)`[2]))
cat("\n--- model WITHOUT first delivery x cesarean interaction: OR (95% profile CI), P ---\n")
out <- data.frame(term = tn$term, OR = round(tn$OR, 3), lo = round(tn$lo, 3), hi = round(tn$hi, 3), p = signif(tn$p.value, 3))
print(out, row.names = FALSE, right = FALSE)
cat("\n--- key terms, full final model vs without the interaction ---\n")
key <- c("Cesarean Section", "1st time pregnant")
for (k in key) { rf <- tf[grepl(paste0("^`?", k), tf$term) & !grepl(":", tf$term), ]; rn <- tn[grepl(paste0("^`?", k), tn$term) & !grepl(":", tn$term), ]
  cat(sprintf("%-20s full: %.2f (%.2f to %.2f)  |  no interaction: %.2f (%.2f to %.2f), P = %.2e\n", k, rf$OR, rf$lo, rf$hi, rn$OR, rn$lo, rn$hi, rn$p.value)) }
write.csv(out, "point2_fig2a_no_firstxCD_interaction.csv", row.names = FALSE)
cat("DONE\n")
