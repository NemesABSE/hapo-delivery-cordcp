# AUTO-GENERATED runner — supplementary forest plots only.
# Sections: 0 Setup, 1 Helpers, 10 Sensitivity, 11 Ancestry, 13 CP95, 14 Parity.
# Bayesian/JointAI sections deliberately excluded (no JAGS needed).
.libPaths(c(normalizePath("~/Rlibs"), .libPaths()))
library(MASS)
library(broom)
library(dplyr)
library(grid)
library(pROC)
library(forestploter)


# rm(list = ls())  # disabled: runner sets libPaths first

# --- Data import (RAW, once) ---
adatok <- read.csv(hapo_file("Hapo_Clinical_FTO_genotype - extended.csv"))
adatok_extra <- read.csv(hapo_file("HAPO_Clinical_with_binary_cordcp.csv"))
adatok$CP_binary <- adatok_extra$Cord_CP_bin[match(adatok$m_mom_geneva_id, adatok_extra$m_mom_geneva_id)]

setwd(normalizePath("figure-relabel/run"))

rm("adatok_extra")
gc()

adatok <- adatok[, -c(1, 2, 4, 5)]

target_variable <- "CP_binary"
id_variable <- "m_mom_geneva_id"

# Regeneration of the Bayesian models (section 9). FALSE BY DEFAULT -> the existing
# .rData is simply loaded. Set it to TRUE ONCE if you want to regenerate (it overwrites
# the imputalt.rData / imputalt_osszevont.rData files), then set it back to FALSE.
REGENERATE_BAYES <- FALSE

# --- Data preparation: factors, GDM, filtering, uniqueness ---
adatok <- adatok %>%
  mutate(
    Ethnicity = factor(Ethnicity, levels = c("White", "Black", "Asian", "Hispanic")),
    m_Hypertens_Binary.None.0.Something.1. = factor(Hypertens.0.no.1.yes.),
    m_firstdelivery = factor(m_OH_PrDel20_1st___1_Other_0__ == 1),
    m_LD_DelType_Binary.Spontan.0.C.Section.1. = factor(ifelse(m_LD_DelType_C_sec__PVN_0_ == 2, 1, 0)),
    b_NN_Gender.Male.1.Female.2. = factor(b_NN_Gender),
    m_Smoker.Yes.1.No.0. = factor(m_Smoker),
    m_FH_DM_Any = factor(m_FH_DM_Any),
    m_FH_HBP_Any = factor(m_FH_HBP_Any),
    m_Drinker = factor(m_Drinker),
    Pre.eclampsia.0.no.1.yes. = factor(Pre.eclampsia.0.no.1.yes.),
    m_Proteinuria = factor(m_Proteinuria),
    m_GDM = ifelse(
      m_FPG_CLC_mmol_l >= 5.1 |
        m_OneHrPG_CLC_mmol_l >= 10 |
        m_TwoHrPG_CLC_mmol_l >= 8.5,
      1, 0
    )
  ) %>%
  filter(
    !is.na(CP_binary),
    !is.na(m_LD_DelType_C_sec__PVN_0_),
    !is.na(Ethnicity),
    !is.na(b_CordCP_ug)
  )

adatok <- adatok %>%
  distinct(.data[[id_variable]], .keep_all = TRUE)

# --- List of explanatory variables (with original, pre-rename names) ---
magyarazo <- c(
  "Ethnicity",
  "m_Hypertens_Binary.None.0.Something.1.",
  "m_firstdelivery",
  "m_LD_DelType_Binary.Spontan.0.C.Section.1.",
  "m_Drinker",
  "m_GDM",
  "b_Gestage",
  "m_HbA1c_percent",
  "b_N_BWMn",
  "b_N_FLMn",
  "b_N_HCMn",
  "b_N_LNGMn",
  "b_N_SSMn",
  "b_N_TRMn",
  "b_NN_Gender.Male.1.Female.2.",
  "m_SBPM_OGTT",
  "m_Smoker.Yes.1.No.0.",
  "m_Age_OGTT",
  "pre_BMI",
  "bmi_at_OGTT",
  "HOMA2_IR",
  "m_FH_DM_Any",
  "m_FH_HBP_Any",
  "WtMgain.till.OGTT",
  "m_Proteinuria",
  "m_HtM_OGTT",
  "Pre.eclampsia.0.no.1.yes.",
  "b_CordPGC_mmol.L"
)

# --- Renaming columns to English labels (data frame + magyarazo list together) ---
adatok <- adatok %>%
  rename(
    `Maternal age at OGTT` = m_Age_OGTT,
    `Maternal GDM` = m_GDM,
    `Age of gestation at delivery` = b_Gestage,
    `Cesarean Section` = m_LD_DelType_Binary.Spontan.0.C.Section.1.,
    Birthweight = b_N_BWMn,
    `Neonatal sex (female)` = b_NN_Gender.Male.1.Female.2.,
    `Maternal BMI at OGTT` = bmi_at_OGTT,
    `1st time pregnant` = m_firstdelivery
  )

# BMI may be raw character -> convert it to numeric (the subgroup analyses work
# directly from `adatok`, and it must be numeric there as well).
adatok$`Maternal BMI at OGTT` <- as.numeric(as.character(adatok$`Maternal BMI at OGTT`))

magyarazo <- recode(
  magyarazo,
  m_Age_OGTT = "Maternal age at OGTT",
  m_GDM = "Maternal GDM",
  b_Gestage = "Age of gestation at delivery",
  m_LD_DelType_Binary.Spontan.0.C.Section.1. = "Cesarean Section",
  b_N_BWMn = "Birthweight",
  b_NN_Gender.Male.1.Female.2. = "Neonatal sex (female)",
  bmi_at_OGTT = "Maternal BMI at OGTT",
  m_firstdelivery = "1st time pregnant"
)

magyarazo <- setdiff(magyarazo, c("m_FCP_ug", "Fasting PG at OGTT"))

# --- "master" data frame: id + magyarazo + target ---
adatok_prepared <- adatok
adatok_master <- adatok_prepared %>%
  select(all_of(c(id_variable, magyarazo, target_variable)))

adatok_master$`Maternal BMI at OGTT` <-
  as.numeric(as.character(adatok_master$`Maternal BMI at OGTT`))

# --- Model variable lists ---
magyarazo_rutinmodel <- c(
  "Maternal age at OGTT",
  "Maternal BMI at OGTT",
  "Ethnicity",
  "1st time pregnant",
  "Cesarean Section",
  "Maternal GDM",
  "Age of gestation at delivery",
  "Birthweight",
  "Neonatal sex (female)"
)
magyarazo_rutinmodel <- intersect(magyarazo_rutinmodel, magyarazo)

magyarazo_bovitett <- magyarazo

# --- CP95 / Parity data (for sections 13 and 14; created here so that they
#     can also be run standalone after Setup) ---
adat_newcp <- read.csv(hapo_file("Hapo_Clinical_FTO_genotype - extended_with_origPrDel20.csv"))
adatok_cp95 <- adatok_master
.idx <- match(adatok_cp95[[id_variable]], adat_newcp$m_mom_geneva_id)
adatok_cp95$CP_binary95 <- adat_newcp$CordCP_95percentile[.idx]
adatok_cp95$m_OH_PrDel20 <- adat_newcp$m_OH_PrDel20[.idx]
adatok_cp95$Parity <- factor(
  ifelse(adatok_cp95$m_OH_PrDel20 >= 3, "3+", as.character(adatok_cp95$m_OH_PrDel20)),
  levels = c("0", "1", "2", "3+")
)



# 1. Helpers -------------------------------------------------------------------

# z-standardization of continuous (non-factor) predictors. Factors are skipped.
standardize_predictors <- function(data, predictors) {
  fakt <- predictors[sapply(data[, predictors, drop = FALSE], is.factor)]
  folyt <- setdiff(predictors, fakt)
  data %>%
    mutate(across(all_of(folyt), ~ as.numeric(scale(.x))))
}

# Building a backtick-quoted formula: every term is wrapped in `...`, joined with "+".
# The `interactions` elements may already contain a `*`/`:` operator with backticks.
build_formula <- function(target, terms, interactions = NULL) {
  rhs <- paste0("`", paste0(terms, collapse = "`+`"), "`")
  if (!is.null(interactions)) {
    rhs <- paste0(rhs, "+", paste0(interactions, collapse = "+"))
  }
  as.formula(paste0(target, "~", rhs))
}

# Safe column renaming: renames only the columns that actually exist.
rename_vars <- function(df, mapping) {
  for (old in names(mapping)) {
    hit <- colnames(df) == old
    if (any(hit)) colnames(df)[hit] <- mapping[[old]]
  }
  df
}

# Same for a character vector (a variable-name list).
rename_vec <- function(x, mapping) {
  for (old in names(mapping)) x[x == old] <- mapping[[old]]
  x
}

# Original -> display names (affects only the existing columns).
lab_map <- c(
  m_Age_OGTT = "Maternal age at OGTT",
  m_GDM = "Maternal GDM",
  b_Gestage = "Age of gestation at delivery",
  m_LD_DelType_Binary.Spontan.0.C.Section.1. = "Cesarean Section",
  b_N_BWMn = "Birthweight",
  b_NN_Gender.Male.1.Female.2. = "Neonatal sex (female)",
  bmi_at_OGTT = "Maternal BMI at OGTT",
  m_firstdelivery = "1st time pregnant",
  m_FCP_nmol_l = "Maternal fasting CP",
  m_FPG_CLC_mmol_l = "Fasting PG at OGTT",
  m_OneHrPG_CLC_mmol_l = "One hour PG at OGTT",
  m_TwoHrPG_CLC_mmol_l = "Two hour PG at OGTT",
  m_OH_PrDel20_faktor = "Parity"
)

# Backtick-quoting interaction term(s) along ":" (for the "joint variable" models).
wrap_terms <- function(terms) {
  sapply(terms, function(x) {
    parts <- strsplit(x, ":", fixed = TRUE)[[1]]
    paste0("`", parts, "`", collapse = ":")
  }, USE.NAMES = FALSE)
}

# p-value formatting (scientific notation).
format_p <- function(p) formatC(p, format = "e", digits = 3)

# Single-model forestplot theme (GLM OR Bayes).
tm_single <- forest_theme(
  base_size = 8,
  ci_Theight = 0.3,
  refline_col = "#63666A",
  ci_col = "#e41a1c",
  footnote_col = "#636363",
  footnote_fontface = "italic",
  legend_name = "Model",
  legend_value = c("GLM", "Bayes"),
  vertline_lty = c("dashed", "dotted"),
  vertline_col = "#bababa",
  row_spacing = unit(4, "mm"),
  core = list(padding = unit(c(4, 3), "mm"))
)

# Two-model forestplot theme (GLM + Bayes together).
tm_double <- forest_theme(
  base_size = 8,
  ci_Theight = 0.3,
  refline_col = "#63666A",
  ci_col = c("#377eb8", "#e41a1c"),
  footnote_col = "#636363",
  footnote_fontface = "italic",
  legend_name = "Model",
  legend_value = c("GLM", "Bayes"),
  vertline_lty = c("dashed", "dotted"),
  vertline_col = c("#d6604d", "#bababa"),
  row_spacing = unit(4, "mm"),
  core = list(padding = unit(c(4, 3), "mm"))
)

# Saving a forestplot to JPG (with the publication parameters).
# Robust device handling: on.exit(dev.off()) guarantees the device is closed even if
# drawing fails -> no "stuck" graphics device is left, which under RStudio would cause
# the "invalid graphics state" (C_palette2) error at the next figure.
save_forest <- function(p, file, width = 3600, height = 3000, res = 500, quality = 100) {
  grDevices::jpeg(file, width = width, height = height, quality = quality, res = res)
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()
  grid::grid.draw(p)
}

# Safe JPG saving of a base-graphics plot (e.g. ROC). on.exit guarantees the dev.off.
save_base_plot <- function(file, expr, width = 2000, height = 2000, res = 500, quality = 100) {
  grDevices::jpeg(file, width = width, height = height, quality = quality, res = res)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(expr)
}

# Single-model OR forestplot: OR (CI) table + plot + JPG.
# est/lower/upper come from the RAW coefficients (the function exponentiates).
or_forest_single <- function(estimate, conf.low, conf.high, p_formatted, labels, file,
                             est_header = "Estimates",
                             or_header = "OR (CI 2.5%-97.5%)",
                             p_header = "P-values",
                             xlim = c(0, 5), ticks = c(0.5, 1, 2, 3),
                             width = 3600, height = 3000,
                             arrow_lab = c("Lower risk", "Higher risk"),
                             fallback_labels = NULL) {
  # The hard-coded labels were made for the earlier variable selection, so it may happen
  # that their length does not match the current model's coefficient count (e.g. b_N_HCMn).
  # In that case we fall back to the raw coefficient names (the figure is still generated).
  if (length(labels) != length(estimate)) {
    warning(sprintf("[%s] label count (%d) != coefficient count (%d); using raw names.",
                    file, length(labels), length(estimate)))
    labels <- if (!is.null(fallback_labels)) fallback_labels else as.character(seq_along(estimate))
  }
  disp <- data.frame(
    Variables = labels,
    est = paste0(round(exp(estimate), 2), " (", round(exp(conf.low), 2),
                 " - ", round(exp(conf.high), 2), ")"),
    or = paste(rep(" ", 40), collapse = " "),
    p = p_formatted,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  names(disp) <- c("Variables", est_header, or_header, p_header)

  # arrow_lab (the Lower/Higher risk arrows) is optional: when NULL we omit it (e.g. for the
  # continuous fold-change figure, where a "risk" label would be misleading).
  forest_args <- list(
    data = disp,
    est = exp(estimate),
    lower = exp(conf.low),
    upper = exp(conf.high),
    ci_column = 3,
    ref_line = 1,
    xlim = xlim,
    ticks_at = ticks,
    theme = tm_single
  )
  if (!is.null(arrow_lab)) forest_args$arrow_lab <- arrow_lab
  p <- do.call(forest, forest_args)
  try(plot(p), silent = TRUE)   # preview; try() guards against a possible RStudio graphics error
  save_forest(p, file, width = width, height = height)
  invisible(p)
}

# "Birth type" combined-interaction model (Cesarean × 1st-time -> 4 categories).
# Same logic for the CP_binary and the CP_binary95 target (only the data/target differ).
birthtype_model <- function(data_master, target) {
  vars <- setdiff(magyarazo_bovitett_r2, "1st time pregnant")
  d <- data_master %>%
    mutate(
      `Birth type` = case_when(
        `Cesarean Section` == 0 & `1st time pregnant` == TRUE ~ "First born, vaginal delivery",
        `Cesarean Section` == 1 & `1st time pregnant` == TRUE ~ "First born, Cesarean section",
        `Cesarean Section` == 0 & `1st time pregnant` == FALSE ~ "Non-first born, vaginal delivery",
        TRUE ~ "Non-first born, Cesarean section"
      ),
      `Birth type` = factor(
        `Birth type`,
        levels = c(
          "First born, vaginal delivery",
          "Non-first born, vaginal delivery",
          "First born, Cesarean section",
          "Non-first born, Cesarean section"
        )
      )
    )
  vars <- c(vars, "Birth type")
  d <- d %>%
    select(all_of(c(id_variable, vars, target))) %>%
    na.omit()
  d <- standardize_predictors(d, vars)

  terms <- c(setdiff(vars, "Cesarean Section"), "Cesarean Section:b_CordPGC_mmol.L")
  terms <- wrap_terms(terms)
  f <- as.formula(paste(target, "~", paste(terms, collapse = " + ")))
  glm(f, family = binomial(), data = d)
}

# --- Forestplot labels (CENTRALIZED, so sections 8/10/11/13 run standalone) ---
# The current models (consistent with the saved selection) contain b_N_HCMn
# (neonatal head circumference) after "Maternal HbA1c at OGTT"; the labels match this.
labels_birthtype <- c(
  "(Intercept)", "Maternal age at OGTT", "Maternal BMI at OGTT", "Ethnicity - Afro-Caribbean",
  "Ethnicity - Thai", "Ethnicity - Hispanic", "Maternal GDM", "Age of gestation at delivery",
  "Birthweight", "Neonatal sex (female)", "Maternal HbA1c at OGTT", "Neonatal head circumference",
  "Neonatal length", "HOMA2-IR at OGTT", "Maternal weight gain (till OGTT)",
  "Neonatal Cord PG concentration", "Non-first born, vaginal delivery", "First born, cesarean delivery",
  "Non-first born, cesarean delivery", "Cesarean delivery * Neonatal Cord PG concentration"
)

# The saved combined Bayesian model (imputalt_osszevont.rData) comes from the OLDER
# selection, so it has 19 terms (without b_N_HCMn) — we use this list for the combined Bayesian figure.
labels_birthtype_bayes <- c(
  "(Intercept)", "Maternal age at OGTT", "Maternal BMI at OGTT", "Ethnicity - Afro-Caribbean",
  "Ethnicity - Thai", "Ethnicity - Hispanic", "Maternal GDM", "Age of gestation at delivery",
  "Birthweight", "Neonatal sex (female)", "Maternal HbA1c at OGTT", "Neonatal length",
  "HOMA2-IR at OGTT", "Maternal weight gain (till OGTT)", "Neonatal Cord PG concentration",
  "Non-first born, vaginal delivery", "First born, cesarean delivery",
  "Non-first born, cesarean delivery", "Cesarean delivery * Neonatal Cord PG concentration"
)

labels_sens1 <- c(
  "(Intercept)", "Maternal age at OGTT", "Maternal BMI at OGTT", "Ethnicity - Afro-Caribbean",
  "Ethnicity - Thai", "Ethnicity - Hispanic", "1st time delivery", "Cesarean delivery",
  "Age of gestation at delivery", "Birthweight", "Neonatal sex (female)",
  "Maternal HbA1c at OGTT", "Neonatal head circumference", "Neonatal length",
  "Maternal weight gain (till OGTT)",
  "Neonatal Cord PG concentration", "Maternal fasting CP", "Fasting PG at OGTT",
  "One hour PG at OGTT", "Two hour PG at OGTT", "1st time delivery * Cesarean delivery",
  "Cesarean delivery * Neonatal Cord PG concentration"
)

labels_sens2 <- c(
  "(Intercept)", "Maternal age at OGTT", "Maternal BMI at OGTT", "Ethnicity - Afro-Caribbean",
  "Ethnicity - Thai", "Ethnicity - Hispanic", "Cesarean delivery", "Maternal GDM",
  "Age of gestation at delivery", "Birthweight", "Neonatal sex (female)",
  "Maternal HbA1c at OGTT", "Neonatal head circumference", "Neonatal length", "HOMA2-IR at OGTT",
  "Maternal weight gain (till OGTT)", "Neonatal Cord PG concentration",
  "Cesarean delivery * Neonatal Cord PG concentration"
)

labels_subpop <- c(
  "(Intercept)", "Maternal age at OGTT", "Maternal BMI at OGTT", "1st time delivery",
  "Cesarean delivery", "Maternal GDM", "Age of gestation at delivery", "Birthweight",
  "Neonatal sex (female)", "Maternal HbA1c at OGTT", "Neonatal head circumference",
  "Neonatal length", "HOMA2-IR at OGTT",
  "Maternal weight gain (till OGTT)", "Neonatal Cord PG concentration",
  "1st time delivery * Cesarean delivery", "Cesarean delivery * Neonatal Cord PG concentration"
)

# Raw coefficient name -> pretty label dictionary (for the new models in section 15, where
# there is no manual label list). pretty_from_terms builds a label list from any model's coef names.
PRETTY <- c(
  "(Intercept)" = "(Intercept)",
  "`Maternal age at OGTT`" = "Maternal age at OGTT",
  "`Maternal BMI at OGTT`" = "Maternal BMI at OGTT",
  "EthnicityBlack" = "Ethnicity - Afro-Caribbean",
  "EthnicityAsian" = "Ethnicity - Thai",
  "EthnicityHispanic" = "Ethnicity - Hispanic",
  "`Maternal GDM`" = "Maternal GDM",
  "`Age of gestation at delivery`" = "Age of gestation at delivery",
  "Birthweight" = "Birthweight",
  "`Neonatal sex (female)`2" = "Neonatal sex (female)",
  "m_HbA1c_percent" = "Maternal HbA1c at OGTT",
  "b_N_HCMn" = "Neonatal head circumference",
  "b_N_LNGMn" = "Neonatal length",
  "b_N_FLMn" = "Neonatal flank skinfold",
  "b_N_SSMn" = "Neonatal subscapular skinfold",
  "b_N_TRMn" = "Neonatal triceps skinfold",
  "HOMA2_IR" = "HOMA2-IR at OGTT",
  "WtMgain.till.OGTT" = "Maternal weight gain (till OGTT)",
  "b_CordPGC_mmol.L" = "Neonatal Cord PG concentration",
  "`1st time pregnant`TRUE" = "1st time delivery",
  "`Cesarean Section`1" = "Cesarean delivery",
  "pre_BMI" = "Maternal pre-pregnancy BMI",
  "`Parity`1" = "Parity: 1 previous delivery",
  "`Parity`2" = "Parity: 2 previous deliveries",
  "`Parity`3+" = "Parity: 3+ previous deliveries",
  # The actual (broom) coef names of the parity model: the Parity main effect without backticks, and
  # the nested Cesarean:Parity terms as "Parity{k}:`Cesarean Section`1" -> cesarean effect per stratum.
  "Parity1" = "Parity: 1 previous delivery",
  "Parity2" = "Parity: 2 previous deliveries",
  "Parity3+" = "Parity: 3+ previous deliveries",
  "Parity0:`Cesarean Section`1" = "Cesarean delivery (Parity 0, firstborn)",
  "Parity1:`Cesarean Section`1" = "Cesarean delivery (Parity 1)",
  "Parity2:`Cesarean Section`1" = "Cesarean delivery (Parity 2)",
  "Parity3+:`Cesarean Section`1" = "Cesarean delivery (Parity 3+)",
  "`1st time pregnant`TRUE:`Cesarean Section`1" = "1st time delivery * Cesarean delivery",
  "`Cesarean Section`1:`1st time pregnant`TRUE" = "1st time delivery * Cesarean delivery",
  "`Cesarean Section`1:b_CordPGC_mmol.L" = "Cesarean delivery * Neonatal Cord PG concentration",
  "b_CordPGC_mmol.L:`Cesarean Section`1" = "Cesarean delivery * Neonatal Cord PG concentration"
)

# Pretty labels from coef names (unknown name -> the raw name is kept).
pretty_from_terms <- function(term_names) {
  out <- PRETTY[term_names]
  out[is.na(out)] <- term_names[is.na(out)]
  unname(out)
}

# --- In case of a missing-data run: na.omit model data frames ---
adatok_bovitett <- na.omit(adatok_master)
adatok_rutin <- na.omit(
  adatok_master[, c(id_variable, magyarazo_rutinmodel, target_variable)]
)

adatok_bovitett <- standardize_predictors(adatok_bovitett, magyarazo_bovitett)
adatok_rutin <- standardize_predictors(adatok_rutin, magyarazo_rutinmodel)

# --- magyarazo_bovitett_r2: the result of the round-1 stepAIC selection. ---
# Section 3 produces it; to run a section standalone we load it from the saved result.
if (!exists("valasztott_valtozok") && file.exists("stepAIC_results.rds")) {
  valasztott_valtozok <- readRDS("stepAIC_results.rds")$selected
}
if (exists("valasztott_valtozok")) {
  magyarazo_bovitett_r2 <- c(magyarazo_rutinmodel, valasztott_valtozok)
}

# --- Shared interaction term list (for the main interaction model and CP95) ---
# Note: it is content-wise IDENTICAL to the original (the term listed twice collapses
# to a single one in the formula anyway, so the model result does not change).
interakcios_tagok_int <- c(
  "`Age of gestation at delivery`*`Birthweight`",
  "`Cesarean Section`*`b_N_LNGMn`",
  "`Cesarean Section`*`WtMgain.till.OGTT`",
  "`Cesarean Section`*`m_HbA1c_percent`",
  "`Ethnicity`*`Maternal GDM`",
  "`Ethnicity`*`b_CordPGC_mmol.L`",
  "`Ethnicity`*`Neonatal sex (female)`",
  "`Ethnicity`*`1st time pregnant`",
  "`Cesarean Section`*`b_CordPGC_mmol.L`",
  "`Cesarean Section`*`1st time pregnant`",
  "`Cesarean Section`*`Maternal age at OGTT`",
  "`Cesarean Section`*`Maternal GDM`",
  "`Ethnicity`*`Maternal age at OGTT`",
  "`Cesarean Section`*`Birthweight`",
  "`Cesarean Section`*`Maternal BMI at OGTT`",
  "`Cesarean Section`*`HOMA2_IR`",
  "`Cesarean Section`*`Ethnicity`",
  "`Ethnicity`*`Maternal BMI at OGTT`"
)

# The two main interactions of the "final model" (for the sensitivity/skinfold/BMI_prepreg/continuous
# models). build_formula backtick-quotes every term; the elements contain an embedded `*`.
final_interactions <- c("Cesarean Section`*`1st time pregnant",
                        "Cesarean Section`*`b_CordPGC_mmol.L")



# 10. Sensitivity --------------------------------------------------------------

# 10a. Glucose and C-peptide values instead of GDM and HOMA
# (labels_sens1 is defined in the Helpers block)
sens_vars <- setdiff(magyarazo_bovitett_r2, c("HOMA2_IR", "Maternal GDM"))
sens_vars <- c(sens_vars, "m_FCP_nmol_l", "m_FPG_CLC_mmol_l",
               "m_OneHrPG_CLC_mmol_l", "m_TwoHrPG_CLC_mmol_l")

uj_valtozok <- adatok %>%
  select(all_of(id_variable), m_FCP_nmol_l, m_FPG_CLC_mmol_l,
         m_OneHrPG_CLC_mmol_l, m_TwoHrPG_CLC_mmol_l)

adatok_sens <- adatok_master %>%
  left_join(uj_valtozok, by = id_variable)
adatok_sens <- rename_vars(adatok_sens, lab_map)
sens_vars <- rename_vec(sens_vars, lab_map)

adatok_sens <- na.omit(adatok_sens[, c(sens_vars, target_variable)])
adatok_sens <- standardize_predictors(adatok_sens, sens_vars)

sens_model_vars <- c(sens_vars,
                     "Cesarean Section`*`1st time pregnant",
                     "Cesarean Section`*`b_CordPGC_mmol.L")
model_sens <- glm(build_formula(target_variable, sens_model_vars),
                  family = "binomial", data = adatok_sens)
summary(model_sens)

co <- broom::tidy(model_sens, conf.int = TRUE)
or_forest_single(
  co$estimate, co$conf.low, co$conf.high, sprintf("%.2e", co$p.value),
  labels_sens1, "sensitivity_1.jpg", fallback_labels = co$term
)

# 10b. Model without birth order
# (labels_sens2 is defined in the Helpers block)
sens2_vars <- setdiff(magyarazo_bovitett_r2, "1st time pregnant")
adatok_sens2 <- na.omit(adatok_master[, c(sens2_vars, target_variable)])
adatok_sens2 <- standardize_predictors(adatok_sens2, sens2_vars)

sens2_model_vars <- c(sens2_vars, "Cesarean Section`*`b_CordPGC_mmol.L")
model_sens2 <- glm(build_formula(target_variable, sens2_model_vars),
                   family = "binomial", data = adatok_sens2)
summary(model_sens2)

co <- broom::tidy(model_sens2, conf.int = TRUE)
or_forest_single(
  co$estimate, co$conf.low, co$conf.high, sprintf("%.2e", co$p.value),
  labels_sens2, "sensitivity_2.jpg", fallback_labels = co$term
)


# 11. Ancestry subgroups -------------------------------------------------------
# White / Black / Asian / Hispanic — each full model + stepAIC + forestplot.
# (labels_subpop is defined in the Helpers block)

fit_subpop <- function(eth, file) {
  d <- na.omit(adatok[adatok$Ethnicity == eth, c(magyarazo_bovitett_r2, target_variable)])
  d <- standardize_predictors(d, magyarazo_bovitett_r2)

  subpop_vars <- c(setdiff(magyarazo_bovitett_r2, "Ethnicity"),
                   "Cesarean Section`*`1st time pregnant",
                   "Cesarean Section`*`b_CordPGC_mmol.L")
  rutin_subpop <- setdiff(magyarazo_rutinmodel, "Ethnicity")

  m_full <- glm(build_formula(target_variable, subpop_vars),
                family = "binomial", data = d)
  m_aic <- stepAIC(
    m_full,
    scope = list(lower = paste0("~`", paste0(rutin_subpop, collapse = "`+`"), "`"))
  )
  print(summary(m_full))
  print(summary(m_aic))
  print(confint(m_aic))
  print(AIC(m_aic))

  co <- broom::tidy(m_full, conf.int = TRUE)
  or_forest_single(
    co$estimate, co$conf.low, co$conf.high, sprintf("%.2e", co$p.value),
    labels_subpop, file, fallback_labels = co$term
  )
  list(data = d, model = m_full, model_aic = m_aic)
}

sp_white <- fit_subpop("White", "sensitivity_eu.jpg")
sp_black <- fit_subpop("Black", "sensitivity_black.jpg")
sp_asian <- fit_subpop("Asian", "sensitivity_asian.jpg")
sp_hisp  <- fit_subpop("Hispanic", "sensitivity_hisp.jpg")

# Compatibility aliases (for the power section, see 12c)
adatok_eu <- sp_white$data;  model_bov_eu <- sp_white$model
adatok_black <- sp_black$data; model_bov_black <- sp_black$model
adatok_hisp <- sp_hisp$data; model_bov_hisp <- sp_hisp$model
model_bov_asian <- sp_asian$model



# 13. Different CP cutoff (CP_binary95) — FIXED ------------------------------
# Requires: Setup + Helpers (adatok_cp95 is built in Setup).
# FIX: the original section rebuilt the already-renamed `adatok` with original
# column names, and CP_binary95 was not included -> it errored out. adatok_cp95
# is now built in Setup, from a fresh CSV read, matched onto adatok_master.
# CAUTION: this section produces a NEW result (it originally did not run).

# 13a. CP95 interaction model (with the logic of the main interaction model)
adatok_cp95_r2 <- na.omit(adatok_cp95[, c(magyarazo_bovitett_r2, "CP_binary95")])
adatok_cp95_r2 <- standardize_predictors(adatok_cp95_r2, magyarazo_bovitett_r2)

model_cp95_int <- glm(
  build_formula("CP_binary95", magyarazo_bovitett_r2, interactions = interakcios_tagok_int),
  family = "binomial", data = adatok_cp95_r2
)
model_cp95_int_aic <- stepAIC(
  model_cp95_int,
  scope = list(lower = paste0("~`", paste0(magyarazo_rutinmodel, collapse = "`+`"),
                              "`+`Ethnicity`*`Maternal BMI at OGTT", "`"))
)
summary(model_cp95_int_aic)
confint(model_cp95_int_aic)
AIC(model_cp95_int_aic)

# 13b. CP95 "Birth type" combined model + forestplot
model_csaszar_cp95 <- birthtype_model(adatok_cp95, "CP_binary95")
summary(model_csaszar_cp95)

co <- broom::tidy(model_csaszar_cp95, conf.int = TRUE)
or_forest_single(
  co$estimate, co$conf.low, co$conf.high, sprintf("%.2e", co$p.value),
  labels_birthtype, "forestplot_csaszar_int_cp95.jpg", fallback_labels = co$term
)

# 13c. Continuous cord C-peptide (log) model — Supp Fig 9b
# The final variables (magyarazo_bovitett_r2) + the two main interactions; outcome: log(b_CordCP_ug).
# linear model -> exp(beta) = per-SD multiplier (fold-change), plotted OR-style.
cp_cont <- adatok %>% select(all_of(id_variable), b_CordCP_ug)
adatok_cont <- adatok_master %>% left_join(cp_cont, by = id_variable)
adatok_cont <- na.omit(adatok_cont[, c(magyarazo_bovitett_r2, "b_CordCP_ug")])
adatok_cont <- adatok_cont[adatok_cont$b_CordCP_ug > 0, ]
adatok_cont <- standardize_predictors(adatok_cont, magyarazo_bovitett_r2)
adatok_cont$logCP <- log(adatok_cont$b_CordCP_ug)

model_cont <- lm(
  build_formula("logCP", c(magyarazo_bovitett_r2, final_interactions)),
  data = adatok_cont
)
summary(model_cont)

co <- broom::tidy(model_cont, conf.int = TRUE)
co <- co[co$term != "(Intercept)" | TRUE, ]  # keep all terms
# The exp(beta)-s of the continuous (log-scale) model are around 1, so a narrow x-axis is needed.
# Larger canvas (5600x3400) so the left labels and right p-values are fully visible; arrow_lab=NULL
# because the "Lower/Higher risk" label would be misleading here (fold-change, not risk).
or_forest_single(
  co$estimate, co$conf.low, co$conf.high, sprintf("%.2e", co$p.value),
  pretty_from_terms(co$term), "supp_fig9_continuous.jpg",
  est_header = "Fold-change per SD (exp beta)", or_header = "Ratio (95% CI)",
  xlim = c(0.5, 1.5), ticks = c(0.5, 0.75, 1, 1.25, 1.5),
  width = 5600, height = 3400, arrow_lab = NULL,
  fallback_labels = co$term
)

# 13d. 95th percentile ROC — Supp Fig 9c
roc_cp95 <- roc(adatok_cp95_r2$CP_binary95, predict(model_cp95_int_aic))
save_base_plot("supp_fig9_roc95.jpg", {
  plot(roc_cp95, col = "#377eb8", lwd = 3, print.auc = TRUE, print.auc.cex = 1.2,
       legacy.axes = TRUE, main = "ROC curve (95th percentile outcome)", grid = TRUE)
  abline(a = 0, b = 1, lty = 2, col = "gray")
})


# 14. Parity — FIXED --------------------------------------------------------
# Requires: Setup + Helpers (adatok_cp95 with the Parity column is built in Setup).
# FIX: the original section ran on corrupted, inherited variable state. Here
# the model is built from adatok_cp95 (with the Parity column), with clean logic.
# DECISION: the target is CP_binary (following the original formula's target_variable2),
# and the filtering is aligned to it. Parity is included instead of "1st time pregnant",
# through the Cesarean Section interactions (Cesarean:Parity, Cesarean:CordPG).
# CAUTION: NEW result (it originally did not run) — please verify the intention.

parity_keep <- setdiff(magyarazo_bovitett_r2, "1st time pregnant")  # Cesarean Section stays in (for the interaction)
adatok_parity <- adatok_cp95 %>%
  select(all_of(c(id_variable, parity_keep, "Parity", "CP_binary"))) %>%
  na.omit()
adatok_parity <- standardize_predictors(adatok_parity, c(parity_keep, "Parity"))

parity_main <- c(setdiff(parity_keep, "Cesarean Section"), "Parity")
parity_int <- c("Cesarean Section`:`Parity", "Cesarean Section`:`b_CordPGC_mmol.L")
model_parity <- glm(
  build_formula("CP_binary", c(parity_main, parity_int)),
  family = "binomial", data = adatok_parity
)
summary(model_parity)

# Supp Fig 8 — parity (categorical 0/1/2/3+) + Cesarean×parity forestplot
co <- broom::tidy(model_parity, conf.int = TRUE)
or_forest_single(
  co$estimate, co$conf.low, co$conf.high, sprintf("%.2e", co$p.value),
  pretty_from_terms(co$term), "supp_fig8_parity.jpg", fallback_labels = co$term
)

