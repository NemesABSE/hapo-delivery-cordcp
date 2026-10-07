## Reviewer 1, Point 2 follow-up (Botond, 2026-09-11): the Fig. 2a Bayesian joint model (JointAI, imputation on the
## 4,951 pairs) refitted WITHOUT the first delivery x cesarean delivery interaction, so that the cesarean coefficient is
## a single birth-order-averaged effect (Bayesian counterpart of point2_fig2a_no_firstxCD_interaction.R, GLM 2.91).
## CD x cord PG and BMI x ancestry interactions kept. Settings as the other rev2 fits: seed 2020, 3 chains, 300
## adaptation, 5,300 iterations, thin 5. Run from this folder.
Sys.setenv(REV2_ROOT = normalizePath("..")); source("../_common/setup.R")
suppressMessages(library(JointAI)); pdf(NULL)
future::plan(future::multisession, workers = 3)
REV <- Sys.getenv("REV2_OUT", unset = ".")
stopifnot(exists("lab_map"), exists("rename_vars"))

d <- rename_vars(adatok_master, lab_map)
d <- d[, c(magyarazo_bovitett_r2, "CP_binary")]
d <- standardize_predictors(d, intersect(magyarazo_bovitett_r2, colnames(d)))
cat("rows:", nrow(d), " (expected 4951) | missing birth order:", sum(is.na(d$`1st time pregnant`)),
    " | missing cesarean:", sum(is.na(d$`Cesarean Section`)), "\n")

ov_terms <- c(magyarazo_bovitett_r2, "Cesarean Section`:`b_CordPGC_mmol.L", "Ethnicity`:`Maternal BMI at OGTT")
f <- paste0("CP_binary~`", paste0(ov_terms, collapse = "`+`"), "`")
clean <- function(x) { x <- gsub(" ", "", x); x <- gsub("_", "", x); x <- gsub("1st", "first", x); gsub("\\(|\\)", "", x) }
f <- clean(f); colnames(d) <- clean(colnames(d))
cat("formula:", f, "\n")

t0 <- Sys.time()
modelimput <- glm_imp(as.formula(f), family = "binomial", data = d, n.adapt = 300, n.iter = 5300, thin = 5,
                      seed = 2020, progress.bar = "none")
cat("fit took", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
save(modelimput, file = file.path(REV, "imputalt_fig2a_no_firstxCD_IMPUTED.rData"))

cl <- modelimput$coef_list$CPbinary; M <- do.call(rbind, lapply(modelimput$MCMC, as.matrix))
B <- M[, cl$coef, drop = FALSE]; colnames(B) <- cl$varname
cat("draws:", nrow(B), " betas:", ncol(B), " data rows:", nrow(modelimput$data), " imputed cells:", sum(is.na(modelimput$data)), "\n")
orq <- function(z) sprintf("%.2f (%.2f to %.2f)  tail-p %.2e", exp(mean(z)), exp(quantile(z, .025)), exp(quantile(z, .975)), 2 * pnorm(-abs(mean(z) / sd(z))))
cat("\n--- all coefficients, posterior mean OR (95% CrI) ---\n")
for (i in seq_len(ncol(B))) cat(sprintf("  %-52s %s\n", colnames(B)[i], orq(B[, i])))
out <- data.frame(varname = colnames(B), OR = exp(colMeans(B)), lo = exp(apply(B, 2, quantile, .025)), hi = exp(apply(B, 2, quantile, .975)),
                  tailp = apply(B, 2, function(z) 2 * pnorm(-abs(mean(z) / sd(z)))))
write.csv(out, "point2_jointai_no_firstxCD_coefs.csv", row.names = FALSE)
cat("\nGelman-Rubin (max Rhat over betas):\n"); print(tryCatch(max(GR_crit(modelimput, multivariate = FALSE)$psrf[, 1], na.rm = TRUE), error = function(e) conditionMessage(e)))
cat("DONE\n")
