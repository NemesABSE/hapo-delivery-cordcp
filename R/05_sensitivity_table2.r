# Sensitivity summary table — how the cesarean delivery and first delivery
# associations behave across every sensitivity and subgroup analysis.
#
# Firneisz's request, 2026-07-21. The difficulty he raised is that some
# analyses re-parameterize one or both of the two variables, so there is no
# single "cesarean OR" to put in a cell. The table is therefore split:
#
#   Block A  both variables keep their binary form -> two directly comparable
#            estimates per analysis, all against the same reference categories
#            (vaginal delivery; non-first delivery)
#   Block B  the analysis re-parameterizes one or both -> one row per level,
#            each naming its own reference category. Block B estimates are NOT
#            comparable with Block A, because the reference differs.
#
# All estimates are complete-case GLM; the Bayesian joint model was only fitted
# for the main model, so there is no Bayesian column to fill for the others.
#
# HOW TO RUN
#   1) Mount the HAPO share so that <HAPO_DATA> resolves.
#   2) Run "0. Setup" and "1. Helpers" of
#      01_main_pipeline.r with the <HAPO_DATA>/
#      paths rewritten to <HAPO_DATA>/ and stepAIC_results.rds in the
#      working directory.
#   3) source() this file. Writes sensitivity_summary_table.csv.
#
# Model definitions are copied verbatim from the sections named in each block
# so the table cannot drift from the figures.

library(MASS); library(broom); library(dplyr)

CD  <- "`Cesarean Section`1"
FTD <- "`1st time pregnant`TRUE"

fmt <- function(co, term) {
  k <- co[co$term == term, ]
  if (!nrow(k)) return(NA_character_)
  sprintf("%.2f (%.2f to %.2f)", exp(k$estimate), exp(k$conf.low), exp(k$conf.high))
}

# Several contrasts stacked into one cell, for the re-parameterized analyses.
fmt_many <- function(co, terms, labels) {
  out <- character(0)
  for (i in seq_along(terms)) {
    v <- fmt(co, terms[i])
    if (!is.na(v)) out <- c(out, paste0(labels[i], " ", v))
  }
  paste(out, collapse = "; ")
}

# One row per scenario. `changed` names the driving assumption or input that
# was varied, which is what makes the table scannable (Firneisz, 2026-07-21).
row_of <- function(model, scenario, changed, item, block,
                   cd_cell = NULL, ftd_cell = NULL, note = "") {
  co <- broom::tidy(model, conf.int = TRUE)
  data.frame(
    Block = block,
    Scenario = scenario,
    `What was varied` = changed,
    n = tryCatch(length(model$y), error = function(e) nobs(model)),
    `Cesarean delivery, OR (95% CI)` =
      if (is.null(cd_cell)) fmt(co, CD) else cd_cell,
    `First delivery, OR (95% CI)` =
      if (is.null(ftd_cell)) fmt(co, FTD) else ftd_cell,
    `Display item` = item,
    Note = note,
    check.names = FALSE, stringsAsFactors = FALSE)
}

rows <- list()

## ---------------------------------------------------------------- Block A
# Final model, complete case (section 4 + 5)
adatok_bovitett_r2 <- na.omit(
  adatok_master[, c(magyarazo_bovitett_r2, target_variable, id_variable)])
adatok_bovitett_r2 <- standardize_predictors(adatok_bovitett_r2,
                                             magyarazo_bovitett_r2)
model_bov_int <- glm(
  build_formula(target_variable, magyarazo_bovitett_r2,
                interactions = interakcios_tagok_int),
  family = "binomial", data = adatok_bovitett_r2)
model_bov_int_aic <- stepAIC(
  model_bov_int, trace = 0,
  scope = list(lower = paste0("~`", paste0(magyarazo_rutinmodel, collapse = "`+`"),
                              "`+`Ethnicity`*`Maternal BMI at OGTT", "`")))
rows[[length(rows) + 1]] <- row_of(model_bov_int_aic,
  "Final model, complete case", "nothing (reference scenario)",
  "Fig. 2a", "A")

# Supplementary Fig. 1 — OGTT glucose and C-peptide instead of GDM and HOMA2-IR
sens_vars <- setdiff(magyarazo_bovitett_r2, c("HOMA2_IR", "Maternal GDM"))
sens_vars <- c(sens_vars, "m_FCP_nmol_l", "m_FPG_CLC_mmol_l",
               "m_OneHrPG_CLC_mmol_l", "m_TwoHrPG_CLC_mmol_l")
uj_valtozok <- adatok %>%
  select(all_of(id_variable), m_FCP_nmol_l, m_FPG_CLC_mmol_l,
         m_OneHrPG_CLC_mmol_l, m_TwoHrPG_CLC_mmol_l)
adatok_sens <- adatok_master %>% left_join(uj_valtozok, by = id_variable)
adatok_sens <- rename_vars(adatok_sens, lab_map)
sens_vars <- rename_vec(sens_vars, lab_map)
adatok_sens <- na.omit(adatok_sens[, c(sens_vars, target_variable)])
adatok_sens <- standardize_predictors(adatok_sens, sens_vars)
model_sens <- glm(build_formula(target_variable,
  c(sens_vars, "Cesarean Section`*`1st time pregnant",
    "Cesarean Section`*`b_CordPGC_mmol.L")),
  family = "binomial", data = adatok_sens)
rows[[length(rows) + 1]] <- row_of(model_sens,
  "Alternative glycemic covariates",
  "GDM and HOMA2-IR replaced by OGTT glucose (0, 60, 120 min) and fasting C-peptide",
  "Supplementary Fig. 1", "A")

# Supplementary Fig. 2 — birth order excluded
sens2_vars <- setdiff(magyarazo_bovitett_r2, "1st time pregnant")
adatok_sens2 <- na.omit(adatok_master[, c(sens2_vars, target_variable)])
adatok_sens2 <- standardize_predictors(adatok_sens2, sens2_vars)
model_sens2 <- glm(build_formula(target_variable,
  c(sens2_vars, "Cesarean Section`*`b_CordPGC_mmol.L")),
  family = "binomial", data = adatok_sens2)
rows[[length(rows) + 1]] <- row_of(model_sens2, "Birth order excluded",
  "first delivery dropped, so cases missing only birth order are retained",
  "Supplementary Fig. 2", "A", ftd_cell = "not in model")

# Supplementary Fig. 5-8 — ancestry strata
fit_sub <- function(eth) {
  d <- na.omit(adatok[adatok$Ethnicity == eth,
                      c(magyarazo_bovitett_r2, target_variable)])
  d <- standardize_predictors(d, magyarazo_bovitett_r2)
  glm(build_formula(target_variable,
        c(setdiff(magyarazo_bovitett_r2, "Ethnicity"),
          "Cesarean Section`*`1st time pregnant",
          "Cesarean Section`*`b_CordPGC_mmol.L")),
      family = "binomial", data = d)
}
for (x in list(c("White", "European", "Supplementary Fig. 5"),
               c("Black", "Afro-Caribbean", "Supplementary Fig. 6"),
               c("Hispanic", "Hispanic", "Supplementary Fig. 7"),
               c("Asian", "Thai", "Supplementary Fig. 8"))) {
  rows[[length(rows) + 1]] <- row_of(fit_sub(x[1]),
    paste("Ancestry stratum:", x[2]),
    "cohort restricted to one ancestry; ancestry terms dropped",
    x[3], "A")
}

# Supplementary Fig. 9c — cord C-peptide as a continuous log-scale outcome
cp_cont <- adatok %>% select(all_of(id_variable), b_CordCP_ug)
adatok_cont <- adatok_master %>% left_join(cp_cont, by = id_variable)
adatok_cont <- na.omit(adatok_cont[, c(magyarazo_bovitett_r2, "b_CordCP_ug")])
adatok_cont <- adatok_cont[adatok_cont$b_CordCP_ug > 0, ]
adatok_cont <- standardize_predictors(adatok_cont, magyarazo_bovitett_r2)
adatok_cont$logCP <- log(adatok_cont$b_CordCP_ug)
model_cont <- lm(build_formula("logCP",
  c(magyarazo_bovitett_r2, final_interactions)), data = adatok_cont)
r <- row_of(model_cont, "Continuous outcome",
  "binary hyperinsulinemia replaced by log cord C-peptide",
  "Supplementary Fig. 9c", "A",
  note = "linear model: exp(beta) is a fold-change, not an odds ratio")
r$n <- nobs(model_cont)
rows[[length(rows) + 1]] <- r

## ---------------------------------------------------------------- Block B
# Fig. 2b — mode of delivery and birth order combined, vaginal firstborn = ref
model_csaszar <- birthtype_model(adatok_master, target_variable)
co_bt <- broom::tidy(model_csaszar, conf.int = TRUE)
bt <- grep("Birth type", co_bt$term, value = TRUE)
bt_lab <- sub(".*`Birth type`", "", bt)
rows[[length(rows) + 1]] <- row_of(model_csaszar,
  "Combined delivery-order variable",
  "mode of delivery and birth order merged into one four-level variable",
  "Fig. 2b", "B",
  cd_cell = fmt_many(co_bt, bt, bt_lab), ftd_cell = "merged into the same variable",
  note = "reference: vaginally delivered firstborn")

# Supplementary Fig. 9a — the same combined variable at the 95th percentile
model_csaszar_cp95 <- birthtype_model(adatok_cp95, "CP_binary95")
co_bt95 <- broom::tidy(model_csaszar_cp95, conf.int = TRUE)
bt95 <- grep("Birth type", co_bt95$term, value = TRUE)
rows[[length(rows) + 1]] <- row_of(model_csaszar_cp95,
  "Combined delivery-order variable, stricter outcome",
  "same merged variable, and the outcome cutoff raised to the 95th percentile",
  "Supplementary Fig. 9a", "B",
  cd_cell = fmt_many(co_bt95, bt95, sub(".*`Birth type`", "", bt95)),
  ftd_cell = "merged into the same variable",
  note = "reference: vaginally delivered firstborn")

# Supplementary Fig. 3 — parity in exact form, with cesarean by parity
parity_keep <- setdiff(magyarazo_bovitett_r2, "1st time pregnant")
adatok_parity <- adatok_cp95 %>%
  select(all_of(c(id_variable, parity_keep, "Parity", "CP_binary"))) %>% na.omit()
adatok_parity <- standardize_predictors(adatok_parity, c(parity_keep, "Parity"))
model_parity <- glm(build_formula("CP_binary",
  c(setdiff(parity_keep, "Cesarean Section"), "Parity",
    "Cesarean Section`:`Parity", "Cesarean Section`:`b_CordPGC_mmol.L")),
  family = "binomial", data = adatok_parity)
co_par <- broom::tidy(model_parity, conf.int = TRUE)
pt_cd <- grep("Parity.*Cesarean", co_par$term, value = TRUE)
pt_vag <- setdiff(grep("^Parity", co_par$term, value = TRUE), pt_cd)
rows[[length(rows) + 1]] <- row_of(model_parity, "Parity in exact form",
  "first delivery (yes/no) replaced by parity 0/1/2/3+, with cesarean by parity",
  "Supplementary Fig. 3", "B",
  cd_cell = fmt_many(co_par, pt_cd, sub("`Cesarean Section`1", "", pt_cd)),
  ftd_cell = fmt_many(co_par, pt_vag, pt_vag),
  note = "reference: vaginal delivery, parity 0")

# Supplementary Fig. 4 — primary versus repeat cesarean delivery
extra <- read.csv(hapo_file("HAPO_Clinical_with_binary_cordcp.csv"))
dt <- extra$m_LD_DelType[match(adatok_cp95[[id_variable]], extra$m_mom_geneva_id)]
adatok_cp95$DelType3 <- factor(dplyr::case_when(
  dt %in% c(1, 2) ~ "Vaginal delivery", dt == 3 ~ "Primary cesarean delivery",
  dt == 4 ~ "Repeat cesarean delivery"),
  levels = c("Vaginal delivery", "Primary cesarean delivery",
             "Repeat cesarean delivery"))
cd_keep <- setdiff(magyarazo_bovitett_r2, "Cesarean Section")
adatok_cd <- adatok_cp95 %>%
  select(all_of(c(id_variable, cd_keep, "DelType3", "CP_binary"))) %>% na.omit()
adatok_cd <- adatok_cd[!(adatok_cd$DelType3 == "Repeat cesarean delivery" &
                         adatok_cd$`1st time pregnant` == TRUE), ]
adatok_cd <- standardize_predictors(adatok_cd, cd_keep)
model_cd3 <- glm(build_formula("CP_binary",
  c(cd_keep, "DelType3", "DelType3`:`b_CordPGC_mmol.L",
    "DelType3`:`1st time pregnant")), family = "binomial", data = adatok_cd)
co_cd3 <- broom::tidy(model_cd3, conf.int = TRUE)
rows[[length(rows) + 1]] <- row_of(model_cd3,
  "Primary versus repeat cesarean delivery",
  "cesarean delivery split into primary and repeat; 10 impossible records excluded",
  "Supplementary Fig. 4", "B",
  cd_cell = fmt_many(co_cd3,
    c("DelType3Primary cesarean delivery", "DelType3Repeat cesarean delivery"),
    c("Primary", "Repeat")),
  note = "reference: vaginal delivery, including instrumental")

out <- do.call(rbind, rows)
print(out, row.names = FALSE)
write.csv(out, "sensitivity_summary_table.csv", row.names = FALSE)
cat("\nrows:", nrow(out), "\n")
