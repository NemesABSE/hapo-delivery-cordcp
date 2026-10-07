## Which generation of the stepwise tables (Supplementary Tables 1 and 2) does the current pipeline reproduce?
## Refits the round-1 and round-2 models from the canonical setup and prints their coefficients.
Sys.setenv(REV2_ROOT = normalizePath("..")); source("../_common/setup.R")
cat("routine:", magyarazo_rutinmodel, sep = " | "); cat("\nselected:", valasztott_valtozok, sep = " | "); cat("\nall candidates:", magyarazo_bovitett, sep = " | "); cat("\n")
tab <- function(m) { co <- broom::tidy(m, conf.int = TRUE); co[, -1] <- round(co[, -1], 3); co$p.value <- sprintf("%.3e", broom::tidy(m)$p.value); print(as.data.frame(co), row.names = FALSE); cat("n =", nobs(m), "\n") }
## round 1: complete cases over every candidate, stepAIC with the routine model as lower bound
d1 <- na.omit(adatok_master[, unique(c(magyarazo_bovitett, target_variable))]); d1 <- standardize_predictors(d1, unique(magyarazo_bovitett))
full <- glm(build_formula(target_variable, unique(magyarazo_bovitett)), family = binomial(), data = d1)
low <- build_formula(target_variable, magyarazo_rutinmodel)
s1 <- MASS::stepAIC(full, scope = list(lower = low), direction = "backward", trace = 0); cat("\n== ROUND 1 ==\n"); tab(s1)
## round 2: complete cases over the round-1 selection, selection repeated
v2 <- all.vars(formula(s1))[-1]; v2 <- gsub("`", "", attr(terms(s1), "term.labels"))
d2 <- na.omit(adatok_master[, c(v2, target_variable)]); d2 <- standardize_predictors(d2, v2)
f2 <- glm(build_formula(target_variable, v2), family = binomial(), data = d2)
s2 <- MASS::stepAIC(f2, scope = list(lower = low), direction = "backward", trace = 0); cat("\n== ROUND 2 ==\n"); tab(s2)
out <- function(m, f) { co <- broom::tidy(m, conf.int = TRUE); write.csv(data.frame(term = co$term, n = nobs(m), co[, -1]), f, row.names = FALSE) }
out(s1, "round1.csv"); out(s2, "round2.csv"); cat("DONE\n")
