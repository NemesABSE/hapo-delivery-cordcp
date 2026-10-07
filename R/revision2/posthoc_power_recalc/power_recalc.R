## Post-hoc power for the cesarean delivery association in the ancestry-stratified fits (Supplementary Table 6, was 9),
## recomputed 2026-09-18. The pipeline (01_main_pipeline.r, power_row) calls pwr.2p.test(h, n = nrow(stratum)), which
## treats the whole stratum as the PER-GROUP size and takes p0 from the whole stratum. Here: (a) the original call
## reproduced, (b) pwr.2p2n.test with the actual cesarean / vaginal group sizes and p0 from the vaginal group,
## (c) the observed power of the adjusted Wald test itself (from beta and SE), which is what the P value implies.
Sys.setenv(REV2_ROOT = normalizePath("..")); source("../_common/setup.R")
library(pwr)
y <- target_variable; lev <- levels(d_cc_s$Ethnicity)
res <- do.call(rbind, lapply(lev, function(l) {
  d <- d_cc %>% filter(Ethnicity == l); vars <- setdiff(magyarazo_bovitett_r2, "Ethnicity")
  ds <- standardize_predictors(d, vars); m <- fit_final(vars, int3[!grepl("Ethnicity", int3)], ds)
  cf <- summary(m)$coefficients["`Cesarean Section`1", ]; beta <- cf["Estimate"]; se <- cf["Std. Error"]; OR <- exp(beta)
  yy <- as.numeric(as.character(d[[y]])) == 1; cd <- as.numeric(as.character(d[["Cesarean Section"]])) == 1
  n <- nrow(d); n_cd <- sum(cd); n_vag <- sum(!cd)
  p_from_or <- function(p0) (OR * p0) / (1 - p0 + OR * p0)
  p0_all <- mean(yy); p0_vag <- mean(yy[!cd])
  pw_orig <- pwr.2p.test(h = ES.h(p_from_or(p0_all), p0_all), n = n, sig.level = 0.05)$power
  h_vag <- ES.h(p_from_or(p0_vag), p0_vag)
  pw_2p2n <- pwr.2p2n.test(h = h_vag, n1 = n_cd, n2 = n_vag, sig.level = 0.05)$power
  z <- abs(beta / se); pw_wald <- pnorm(z - qnorm(0.975)) + pnorm(-z - qnorm(0.975))
  data.frame(ancestry = l, N = n, n_cesarean = n_cd, n_vaginal = n_vag, events = sum(yy), events_cesarean = sum(yy & cd),
             events_vaginal = sum(yy & !cd), prev_all = p0_all, prev_vaginal = p0_vag, prev_cesarean = mean(yy[cd]),
             beta = beta, SE = se, OR = OR, p_wald = cf["Pr(>|z|)"], h_original = ES.h(p_from_or(p0_all), p0_all),
             power_original_pwr2p_nN = pw_orig, h_vaginal_ref = h_vag, power_2p2n_actual_groups = pw_2p2n,
             power_observed_wald = pw_wald, row.names = NULL) }))
options(width = 250); print(res, digits = 3, row.names = FALSE)
write.csv(res, "power_recalc.csv", row.names = FALSE)
cat("\nDONE\n")
