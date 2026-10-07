## Reviewer 1, Point 3: birth-order missingness. (1) stratum of Supp Table 5, (2) imputation-model
## structure of the JointAI joint model, (3) delta-adjusted tipping-point (MNAR) sensitivity,
## (4) Thai complete-case subset vs excluded subset.
## Run from this folder: Rscript point3.R > point3.log 2>&1
Sys.setenv(REV2_ROOT = normalizePath("..")); source("../_common/setup.R")
set.seed(2020)

## ---- (1) Supplementary Table 5 refit, Thai stratum ---------------------------------------------
cat("\n[3.1] Supplementary Table 5: logistic model of birth-order missingness, Thai stratum\n")
thai <- adatok_master %>% filter(Ethnicity == "Asian")
thai$miss_bo <- as.integer(is.na(thai$`1st time pregnant`))
cat("Thai n =", nrow(thai), " missing birth order =", sum(thai$miss_bo), "(", round(100 * mean(thai$miss_bo), 1), "% )\n")
m5 <- glm(miss_bo ~ CP_binary + `Maternal age at OGTT` + `Maternal BMI at OGTT` + I(Birthweight / 100), family = binomial(), data = thai)
cat("n in model:", nobs(m5), "\n"); print(orci(m5), digits = 3, row.names = FALSE)
cat("(published: hyperinsulinemia OR 2.54, 1.65 to 3.95; age 0.87; BMI 0.96; birthweight/100 g 0.97)\n")
cat("\n[3.4] hyperinsulinemia frequency, Thai: missing vs known birth order\n")
print(with(thai, tapply(CP_binary, miss_bo, function(x) sprintf("%.1f%% (n = %d)", 100 * mean(x, na.rm = TRUE), sum(!is.na(x))))))
thai_cc <- d_cc %>% filter(Ethnicity == "Asian")
cat("Thai complete-case set of the final model: n =", nrow(thai_cc), " hyperinsulinemia =", round(100 * mean(thai_cc$CP_binary), 1), "%\n")

## ---- (2) structure of the joint model (what imputes birth order) ---------------------------------
cat("\n[3.2] JointAI joint model (imputalt.rData): covariate sub-models\n")
e <- new.env(); load(file.path(hapo_bayes_dir(), "imputalt.rData"), envir = e); m9 <- get(ls(e)[1], envir = e)
cat("analysis model:", paste(m9$coef_list$CPbinary$varname[-1], collapse = " + "), "\n")
for (nm in setdiff(names(m9$coef_list), "CPbinary")) cat("  ", nm, "~", paste(m9$coef_list[[nm]]$varname[-1], collapse = " + "), "\n")
na <- colSums(is.na(m9$data)); cat("missing cells per variable:\n"); print(na[na > 0])
cat("Note: in a joint model the missing covariate values are drawn from their full conditional,\n",
    "which is the product of the covariate sub-model and the OUTCOME likelihood, so the imputation\n",
    "of birth order is conditioned on hyperinsulinemia even though the outcome is not a regressor\n",
    "in the covariate sub-model (Erler et al. 2016 / JointAI documentation).\n")

## ---- (3) delta-adjusted tipping-point sensitivity ------------------------------------------------
cat("\n[3.3] Delta-adjusted (pattern-mixture) tipping-point analysis for missing birth order\n")
oth <- setdiff(magyarazo_bovitett_r2, "1st time pregnant")
d_tp <- adatok_master %>% select(all_of(c(id_variable, oth, "1st time pregnant", target_variable))) %>%
  filter(if_all(all_of(c(oth, target_variable)), ~ !is.na(.x)))
cat("rows complete for every final-model variable except birth order:", nrow(d_tp), " of which birth order missing:", sum(is.na(d_tp$`1st time pregnant`)), "\n")
print(table(ancestry = d_tp$Ethnicity, missing_birth_order = is.na(d_tp$`1st time pregnant`)))
d_tp_s <- standardize_predictors(d_tp, oth)
d_tp_s$first01 <- as.integer(as.logical(d_tp_s$`1st time pregnant`))
obs <- !is.na(d_tp_s$first01)
f_imp <- as.formula(paste("first01 ~", paste0("`", oth, "`", collapse = " + "), "+", target_variable, "+ `Cesarean Section`:b_CordPGC_mmol.L"))
m_imp <- glm(f_imp, family = binomial(), data = d_tp_s[obs, ])
cat("\nimputation model for birth order (fitted on the", sum(obs), "observed rows): predictors = all final-model covariates + hyperinsulinemia + CD x cord PG\n")
oi <- orci(m_imp); print(oi[grepl("CP_binary|Cesarean|Ethnicity|age", oi$term), ], digits = 3, row.names = FALSE)
tmp <- d_tp_s; tmp$first01[!obs] <- 0L
X_mis <- model.matrix(f_imp, data = tmp[!obs, ])
bhat <- coef(m_imp); Vhat <- vcov(m_imp)
rmvn <- function(mu, S) as.numeric(mu + t(chol(S)) %*% rnorm(length(mu)))
pool <- function(est, var) { m <- length(est); q <- mean(est); U <- mean(var); Bv <- var(est); Tt <- U + (1 + 1/m) * Bv
  df <- (m - 1) * (1 + U / ((1 + 1/m) * Bv))^2; c(q, sqrt(Tt), df) }
deltas <- c(-3, -2, -1, -0.5, 0, 0.5, 1, 2, 3); M <- 25
f_an <- build_formula(target_variable, c(magyarazo_bovitett_r2, int3))
se_of <- function(cm, i) (log(cm[i, 3]) - log(cm[i, 2])) / 3.92
rows <- list()
for (dl in deltas) {
  E <- list(); pf <- numeric(M)
  for (k in 1:M) {
    pi_k <- plogis(as.numeric(X_mis %*% rmvn(bhat, Vhat)) + dl)
    dk <- d_tp_s; dk$first01[!obs] <- rbinom(sum(!obs), 1, pi_k); pf[k] <- mean(dk$first01[!obs])
    dk$`1st time pregnant` <- factor(dk$first01 == 1)
    cm <- birthtype_contrasts(glm(f_an, family = binomial(), data = dk))
    E[[k]] <- c(first = -log(cm[1, 1]), se_first = se_of(cm, 1), cd_nf = log(cm[4, 1]), se_cd_nf = se_of(cm, 4),
                cd_f = log(cm[5, 1]), se_cd_f = se_of(cm, 5), nf_cd = log(cm[3, 1]), se_nf_cd = se_of(cm, 3))
  }
  E <- do.call(rbind, E)
  P <- function(a, s) { r <- pool(E[, a], E[, s]^2); fmt(c(exp(r[1]), exp(r[1] - qt(.975, r[3]) * r[2]), exp(r[1] + qt(.975, r[3]) * r[2]))) }
  rows[[length(rows) + 1]] <- data.frame(delta = dl, imputed_firstborn_pct = round(100 * mean(pf), 1),
    OR_first_delivery = P("first", "se_first"), OR_CD_within_nonfirst = P("cd_nf", "se_cd_nf"),
    OR_CD_within_first = P("cd_f", "se_cd_f"), OR_nonfirst_CD_vs_first_vaginal = P("nf_cd", "se_nf_cd"))
}
tp <- do.call(rbind, rows)
cat("\nobserved firstborn % among known birth order:", round(100 * mean(d_tp_s$first01[obs]), 1), "\n")
cat("delta = log-odds shift added to the MAR imputation model for P(firstborn) in the", sum(!obs), "pairs with missing birth order;",
    "M =", M, "proper imputations (parameters drawn from the asymptotic posterior), Rubin's rules; delta = 0 is MAR.\n\n")
print(tp, right = FALSE, row.names = FALSE)
write.csv(tp, "point3_tipping_point.csv", row.names = FALSE)
cat("\nDONE\n")
