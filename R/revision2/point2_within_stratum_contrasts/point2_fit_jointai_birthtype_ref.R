## Reviewer 1, Point 2: the Fig. 3 Bayesian joint model (four-level Birth type + cesarean x cord glucose +
## BMI x ancestry, imputed on the 4,951 pairs) REFITTED with the reference level of Birth type set to
## "Non-first born, vaginal delivery", so that the cesarean vs vaginal odds ratio WITHIN NON-FIRSTBORNS is a
## directly estimated coefficient (the "Non-first born, Cesarean section" row), not a derived contrast.
## Identical settings to the option-B fit (seed 2020, 3 chains, 300 adaptation, 5,300 iterations, thin 5).
## Reference level via env REV2_REF (default below); output file name carries the reference.
## Run from this folder: Rscript point2_fit_jointai_birthtype_ref.R > point2_fit_jointai_birthtype_ref.log 2>&1
Sys.setenv(REV2_ROOT = normalizePath("..")); source("../_common/setup.R")
suppressMessages(library(JointAI)); pdf(NULL)
future::plan(future::multisession, workers = 3)     # 3 chains on 3 cores (JointAI >= 1.0 parallelises via future)
REV <- Sys.getenv("REV2_OUT", unset = ".")
REF <- Sys.getenv("REV2_REF", unset = "Non-first born, vaginal delivery")
TAG <- gsub("[^A-Za-z]", "", REF)
stopifnot(exists("lab_map"), exists("rename_vars"))

d <- rename_vars(adatok_master, lab_map)
d <- d[, c(magyarazo_bovitett_r2, "CP_binary")]
d <- standardize_predictors(d, intersect(magyarazo_bovitett_r2, colnames(d)))
cat("rows:", nrow(d), " (expected 4951) | missing birth order:", sum(is.na(d$`1st time pregnant`)),
    " | missing cesarean:", sum(is.na(d$`Cesarean Section`)), "\n")
d$`Birth type` <- factor(dplyr::case_when(
  d$`Cesarean Section` == 0 & d$`1st time pregnant` == TRUE  ~ "First born, vaginal delivery",
  d$`Cesarean Section` == 1 & d$`1st time pregnant` == TRUE  ~ "First born, Cesarean section",
  d$`Cesarean Section` == 0 & d$`1st time pregnant` == FALSE ~ "Non-first born, vaginal delivery",
  d$`Cesarean Section` == 1 & d$`1st time pregnant` == FALSE ~ "Non-first born, Cesarean section",
  TRUE ~ NA_character_),                                    # NA birth order stays NA (no catch-all)
  levels = c("First born, vaginal delivery", "Non-first born, vaginal delivery",
             "First born, Cesarean section", "Non-first born, Cesarean section"))
d$`Birth type` <- relevel(d$`Birth type`, ref = REF)
cat("Birth type (reference =", REF, "):\n"); print(table(d$`Birth type`, useNA = "ifany"))
d$`1st time pregnant` <- NULL                                # replaced by Birth type; CD kept for the interaction only

ov_terms <- c(setdiff(magyarazo_bovitett_r2, c("Cesarean Section", "1st time pregnant")), "Birth type",
              "Cesarean Section`:`b_CordPGC_mmol.L", "Ethnicity`:`Maternal BMI at OGTT")
f <- paste0("CP_binary~`", paste0(ov_terms, collapse = "`+`"), "`")
clean <- function(x) { x <- gsub(" ", "", x); x <- gsub("_", "", x); x <- gsub("1st", "first", x); gsub("\\(|\\)", "", x) }
f <- clean(f); colnames(d) <- clean(colnames(d)); levels(d$Birthtype) <- gsub(" ", "", levels(d$Birthtype))
cat("formula:", f, "\n")

t0 <- Sys.time()
modelimput <- glm_imp(as.formula(f), family = "binomial", data = d, n.adapt = 300, n.iter = 5300, thin = 5,
                      seed = 2020, progress.bar = "none")
cat("fit took", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
save(modelimput, file = file.path(REV, paste0("imputalt_birthtype_withBMIxAnc_ref", TAG, "_IMPUTED.rData")))

cl <- modelimput$coef_list$CPbinary; M <- do.call(rbind, lapply(modelimput$MCMC, as.matrix))
B <- M[, cl$coef, drop = FALSE]; colnames(B) <- cl$varname
cat("draws:", nrow(B), " betas:", ncol(B), " data rows:", nrow(modelimput$data), " imputed cells:", sum(is.na(modelimput$data)), "\n")
orq <- function(z) sprintf("%.2f (%.2f to %.2f)  tail-p %.2e", exp(mean(z)), exp(quantile(z, .025)), exp(quantile(z, .975)), 2 * pnorm(-abs(mean(z) / sd(z))))
cat("\n--- all coefficients, posterior mean OR (95% CrI) ---\n")
for (i in seq_len(ncol(B))) cat(sprintf("  %-52s %s\n", colnames(B)[i], orq(B[, i])))
out <- data.frame(varname = colnames(B), OR = exp(colMeans(B)), lo = exp(apply(B, 2, quantile, .025)), hi = exp(apply(B, 2, quantile, .975)),
                  tailp = apply(B, 2, function(z) 2 * pnorm(-abs(mean(z) / sd(z)))))
write.csv(out, paste0("point2_jointai_birthtype_ref", TAG, "_coefs.csv"), row.names = FALSE)
cat("\nGelman-Rubin (max Rhat over betas):\n"); print(tryCatch(max(GR_crit(modelimput, multivariate = FALSE)$psrf[, 1], na.rm = TRUE), error = function(e) conditionMessage(e)))
cat("DONE\n")
