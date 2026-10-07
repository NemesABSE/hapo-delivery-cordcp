## Supp Table 6 rework (2026-10-05): observed Wald power of the GLOBAL cesarean delivery coefficient (model without the
## cesarean delivery x first delivery interaction, as in Table 2) in the full complete-case sample and in each ancestry
## stratum. Same fits as table2_global_cd_0924/table2_global_cd.R sections 1 and 5.
## Rscript power_global.R > power_global.log 2>&1
Sys.setenv(REV2_ROOT = normalizePath("..")); source("../_common/setup.R")
drop_first_cd <- function(ints) ints[!grepl("1st time pregnant", ints)]
b_cd <- "`Cesarean Section`1"
pw <- function(b, se, alpha = 0.05) { z <- abs(b / se); cr <- qnorm(1 - alpha / 2); pnorm(z - cr) + pnorm(-z - cr) }
row_of <- function(lab, m, d) { s <- summary(m)$coefficients[b_cd, ]; y <- d[[target_variable]]; cd <- d[["Cesarean Section"]] == 1
  ev <- as.numeric(as.character(y)) == 1 | y == TRUE
  data.frame(ancestry = lab, N = nobs(m), n_cesarean = sum(cd), n_vaginal = sum(!cd), events = sum(ev),
    beta = s[1], SE = s[2], OR = exp(s[1]), p_wald = s[4], power_observed_wald = pw(s[1], s[2]), row.names = NULL) }
res <- list(row_of("All", fit_final(ints = drop_first_cd(int3)), d_cc_s))
for (l in c("White", "Black", "Hispanic", "Asian")) {
  d <- na.omit(adatok[adatok$Ethnicity == l, c(magyarazo_bovitett_r2, target_variable)])
  vars <- setdiff(magyarazo_bovitett_r2, "Ethnicity"); d <- standardize_predictors(d, vars)
  mg <- glm(build_formula(target_variable, c(vars, drop_first_cd(final_interactions))), family = binomial(), data = d)
  res[[length(res) + 1]] <- row_of(l, mg, d)
}
out <- do.call(rbind, res); print(out, digits = 4, row.names = FALSE)
write.csv(out, "power_global.csv", row.names = FALSE); cat("DONE\n")
