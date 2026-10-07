## Supplementary Table 8 (birth-order missingness, Thai stratum) refitted on the current data (maternal BMI at OGTT).
## Same model as pipeline section 12b and point3_missing_birth_order/point3.R; also the version with pre-pregnancy BMI,
## to show where the printed 2.54 came from.
Sys.setenv(REV2_ROOT = normalizePath("..")); source("../_common/setup.R")
thai <- adatok_master %>% filter(Ethnicity == "Asian"); thai$miss_bo <- as.integer(is.na(thai$`1st time pregnant`))
cat("Thai n =", nrow(thai), " missing birth order =", sum(thai$miss_bo), "\n")
m <- glm(miss_bo ~ CP_binary + `Maternal age at OGTT` + `Maternal BMI at OGTT` + I(Birthweight / 100), family = binomial(), data = thai)
o <- orci(m); o$n <- nobs(m); print(o, digits = 6, row.names = FALSE); write.csv(o, "table8_refit.csv", row.names = FALSE)
cat("BMI-like columns:", grep("BMI|bmi", names(adatok), value = TRUE), "\n")
## which BMI gives the printed 2.54 (1.65 to 3.95)?
raw <- adatok_prepared %>% filter(Ethnicity == "Asian"); thai$pre_BMI <- raw$pre_BMI[match(thai[[id_variable]], raw[[id_variable]])]
m0 <- glm(miss_bo ~ CP_binary + `Maternal age at OGTT` + pre_BMI + I(Birthweight / 100), family = binomial(), data = thai)
cat("\npre-pregnancy BMI version, n =", nobs(m0), "\n"); print(orci(m0), digits = 4, row.names = FALSE)
