## Supp Table 6, 2026-10-06: observed Wald power for the three cesarean delivery contrasts in the full complete-case sample
## and in each ancestry stratum: global (model without cesarean delivery x first delivery), and cesarean versus vaginal
## delivery within non-firstborns and within firstborns (interaction model; firstborn = main effect + interaction, delta method).
## Rscript power_three_contrasts.R > power_three_contrasts.log 2>&1
Sys.setenv(REV2_ROOT = normalizePath("..")); source("../_common/setup.R")
drop_first_cd <- function(ints) ints[!grepl("1st time pregnant", ints)]
b_cd <- "`Cesarean Section`1"
pw <- function(z, alpha = 0.05) { z <- abs(z); cr <- qnorm(1 - alpha / 2); pnorm(z - cr) + pnorm(-z - cr) }
one <- function(lab, contrast, n, b, se) data.frame(ancestry = lab, N = n, contrast = contrast, beta = b, SE = se, OR = exp(b),
  p_wald = 2 * pnorm(-abs(b / se)), power_observed_wald = pw(b / se))
rows_of <- function(lab, m, mg) {
  bi <- find_term(m, "Cesarean Section`1:`1st time pregnant`TRUE$|1st time pregnant`TRUE:`Cesarean Section`1$"); V <- vcov(m); b <- coef(m)
  rbind(one(lab, "Global", nobs(mg), coef(mg)[b_cd], sqrt(vcov(mg)[b_cd, b_cd])),
        one(lab, "Non-firstborns", nobs(m), b[b_cd], sqrt(V[b_cd, b_cd])),
        one(lab, "Firstborns", nobs(m), b[b_cd] + b[bi], sqrt(V[b_cd, b_cd] + V[bi, bi] + 2 * V[b_cd, bi])))
}
res <- list(rows_of("All", fit_final(), fit_final(ints = drop_first_cd(int3))))
for (l in c("White", "Black", "Hispanic", "Asian")) {
  d <- na.omit(adatok[adatok$Ethnicity == l, c(magyarazo_bovitett_r2, target_variable)])
  vars <- setdiff(magyarazo_bovitett_r2, "Ethnicity"); d <- standardize_predictors(d, vars)
  m <- glm(build_formula(target_variable, c(vars, final_interactions)), family = binomial(), data = d)
  mg <- glm(build_formula(target_variable, c(vars, drop_first_cd(final_interactions))), family = binomial(), data = d)
  res[[length(res) + 1]] <- rows_of(l, m, mg)
}
out <- do.call(rbind, res); rownames(out) <- NULL; print(out, digits = 4); write.csv(out, "power_three_contrasts.csv", row.names = FALSE); cat("DONE\n")
