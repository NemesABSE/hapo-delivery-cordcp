## Which BMI reproduces the published Supp Table 5 (hyperinsulinemia OR 2.54, 1.65 to 3.95)?
Sys.setenv(REV2_ROOT = normalizePath("..")); source("../_common/setup.R")
thai <- adatok_master %>% filter(Ethnicity == "Asian"); thai$miss_bo <- as.integer(is.na(thai$`1st time pregnant`))
raw <- adatok_prepared %>% filter(Ethnicity == "Asian")
thai$pre_BMI <- raw$pre_BMI[match(thai[[id_variable]], raw[[id_variable]])]
for (v in c("`Maternal BMI at OGTT`", "pre_BMI")) {
  m <- glm(as.formula(paste("miss_bo ~ CP_binary + `Maternal age at OGTT` +", v, "+ I(Birthweight/100)")), family = binomial(), data = thai)
  cat("\nBMI variable:", v, " n =", nobs(m), "\n"); print(orci(m), digits = 3, row.names = FALSE)
}
