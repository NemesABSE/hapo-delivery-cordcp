## Point 4: counts of m_LD_DelType (1 spontaneous vaginal, 2 instrumental, 3 primary CS, 4 repeat CS) in the analysis
## set (n = 4,951) and in the complete-case set of the final model (n = 3,482).
Sys.setenv(REV2_ROOT = normalizePath("..")); source("../_common/setup.R")
cat("\ncolumns matching DelType in adatok_master:", grep("DelType", names(adatok_master), value = TRUE), "\n")
raw <- read.csv(hapo_file("HAPO_Clinical_with_binary_cordcp.csv"), check.names = FALSE)
cat("HAPO_Clinical_with_binary_cordcp.csv: rows", nrow(raw), " DelType cols:", grep("DelType", names(raw), value = TRUE), " id col present:", id_variable %in% names(raw), "\n")
cat("duplicate ids in the clinical file:", sum(duplicated(raw[[id_variable]])), "\n")
.i <- match(adatok_master[[id_variable]], raw[[id_variable]])   # same first-match convention as supp_fig4_primary_repeat_CD.r
cat("matched:", sum(!is.na(.i)), "of", nrow(adatok_master), "\n")
src <- data.frame(adatok_master[, id_variable, drop = FALSE], m_LD_DelType = raw$m_LD_DelType[.i], check.names = FALSE)
v <- grep("DelType", names(src), value = TRUE)[1]
cat("\nAnalysis set (n =", nrow(adatok_master), "):\n"); print(table(src[[v]], useNA = "ifany"))
cc_ids <- d_cc[[id_variable]]
cat("\nComplete-case set (n =", length(cc_ids), "):\n"); print(table(src[[v]][src[[id_variable]] %in% cc_ids], useNA = "ifany"))
cat("\nby ancestry, complete-case:\n"); print(table(adatok_master$Ethnicity[match(cc_ids, adatok_master[[id_variable]])], src[[v]][match(cc_ids, src[[id_variable]])]))
cat("\nhyperinsulinemia by DelType, complete-case:\n"); print(table(src[[v]][match(cc_ids, src[[id_variable]])], d_cc[[target_variable]]))
cat("DONE\n")
