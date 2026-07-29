# ==============================================================================
#  Cord CP — full ethnicity, cutoff analysis  (COMPACT, REFACTORED VERSION)
# ------------------------------------------------------------------------------
#  Original: "Cord_CP_full_ethinc_cutoff_v260629 _ review.r"
#  Goal: more compact, error-free code; the sections can be run INDIVIDUALLY after
#        running the "0. Setup" + "1. Helpers" block.
#  IMPORTANT: the numerical RESULTS of the working sections are unchanged. The detailed
#          refactoring log is in the steps.txt file.
#
#  Running:
#    1) Run the "0. Setup" and "1. Helpers" block.
#    2) After that most sections can be run STANDALONE (after Setup), because every
#       key constant (magyarazo_bovitett_r2, labels_*, interakcios_tagok_int,
#       adatok_cp95) is already defined in the Setup/Helpers block.
#
#  DEPENDENCY CHAIN (what can be run after what):
#    Setup + Helpers  -> MANDATORY for everything else.
#    2, 3, 8, 10, 11, 13, 14, 15  -> only Setup+Helpers is needed for them.
#    4  -> Setup (magyarazo_bovitett_r2). Produces: adatok_bovitett_r2.
#    5  -> 4 (adatok_bovitett_r2). Produces: model_bov_int_aic, ROC.jpg.
#    6  -> 4 (adatok_bovitett_r2) + imputalt.rData. Loads: modelimput.
#    7  -> 5 (model_bov_int_aic) + 6 (modelimput).
#    12 -> the power sub-block (12c) REQUIRES section 11 (ancestry models: sp_*, model_bov_*).
#          12a/12b (Thai distribution, missingness) only require Setup.
#    9  -> REGENERATION of the Bayesian models. Runs only if REGENERATE_BAYES == TRUE.
# ==============================================================================


# 0. Setup ---------------------------------------------------------------------

library(tableone)
library(MASS)
library(performance)
library(DHARMa)
# install.packages('JointAI')
library(JointAI)
library(sensitivity)
library(broom)
library(dplyr)
library(grid)
library(ROCR)
library(pROC)
library(forestploter)
library(flextable)
library(officer)
library(pwr)

rm(list = ls())

# --- Data import (RAW, once) ---
adatok <- read.csv(hapo_file("Hapo_Clinical_FTO_genotype - extended.csv"))
adatok_extra <- read.csv(hapo_file("HAPO_Clinical_with_binary_cordcp.csv"))
adatok$CP_binary <- adatok_extra$Cord_CP_bin[match(adatok$m_mom_geneva_id, adatok_extra$m_mom_geneva_id)]

setwd(Sys.getenv("HAPO_OUT", unset = "."))

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
  # Birth order is deliberately absent from `vars`, so a missing value in it cannot
  # be caught by the na.omit() below; such records used to fall through the case_when
  # catch-all and be labelled "Non-first born, Cesarean section" whatever their actual
  # mode of delivery. Drop them up front so this model is complete-case (n = 3,482)
  # like every other model in the paper.
  d <- data_master %>%
    filter(!is.na(`1st time pregnant`)) %>%
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


# 2. Descriptive statistics ----------------------------------------------------

sapply(adatok_master[, magyarazo], function(x)
  if (is.factor(x) | is.character(x)) length(unique(x)) else NA)

leiro <- CreateTableOne(
  vars = magyarazo,
  strata = target_variable,
  data = adatok_master
)

leiro2 <- CreateTableOne(
  vars = magyarazo_rutinmodel,
  data = adatok_master
)


# 3. Modelling — round 1 -------------------------------------------------------

# Routine model
model_rutin <- glm(build_formula(target_variable, magyarazo_rutinmodel),
                   family = "binomial", data = adatok_rutin)
model_rutin_aic <- stepAIC(model_rutin)

summary(model_rutin)
summary(model_rutin_aic)

# Extended model (lower scope = routine model)
model_bov <- glm(build_formula(target_variable, magyarazo_bovitett),
                 family = "binomial", data = adatok_bovitett)
model_bov_aic <- stepAIC(
  model_bov,
  scope = list(lower = paste0("~`", paste0(magyarazo_rutinmodel, collapse = "`+`"), "`"))
)

summary(model_bov)
summary(model_bov_aic)
confint(model_bov_aic)
AIC(model_bov_aic)

# Selected / dropped variables
valasztott_valtozok <- attr(terms(model_bov_aic), "term.labels")
valasztott_valtozok <- gsub("`", "", valasztott_valtozok)
valasztott_valtozok

kiesett_valtozok <- setdiff(magyarazo_bovitett, valasztott_valtozok)
kiesett_valtozok

saveRDS(
  list(selected = valasztott_valtozok, removed = kiesett_valtozok),
  "stepAIC_results.rds"
)

# Variable list used for sections 4–14 (from the round-1 result)
magyarazo_bovitett_r2 <- c(magyarazo_rutinmodel, valasztott_valtozok)

# Results table (AIC_round1.csv — saving intentionally commented out, as originally)
coefs <- summary(model_bov_aic)$coefficients
confint_vals <- confint(model_bov_aic)
results_df <- data.frame(
  Variable = rownames(coefs),
  Estimate = round(coefs[, "Estimate"], 3),
  Std.Error = round(coefs[, "Std. Error"], 3),
  z.value = round(coefs[, "z value"], 3),
  p.value = sapply(coefs[, "Pr(>|z|)"], format_p),
  CI.Lower = round(confint_vals[, 1], 3),
  CI.Upper = round(confint_vals[, 2], 3),
  stringsAsFactors = FALSE
)
# write.csv2(results_df, "AIC_round1.csv", row.names = FALSE)

# Routine model on the patient set used in the extended one (id-based matching)
adatok3 <- adatok_rutin[adatok_rutin[[id_variable]] %in% adatok_bovitett[[id_variable]], ]
model_rutin_sub <- glm(build_formula(target_variable, magyarazo_rutinmodel),
                       family = "binomial", data = adatok3)
model_rutin_sub_aic <- stepAIC(model_rutin_sub)

summary(model_rutin_sub)
confint(model_rutin_sub)
summary(model_rutin_sub_aic)


# 4. Round 2 — variable selection ----------------------------------------------
# Requires: Setup + Helpers (magyarazo_bovitett_r2). Produces: adatok_bovitett_r2.

# Version without FCP and FPG: routine + round-1 selected variables
adatok_bovitett_r2 <- na.omit(
  adatok_master[, c(magyarazo_bovitett_r2, target_variable, id_variable)]
)
adatok_bovitett_r2 <- standardize_predictors(adatok_bovitett_r2, magyarazo_bovitett_r2)

model_bov_r2 <- glm(build_formula(target_variable, magyarazo_bovitett_r2),
                    family = "binomial", data = adatok_bovitett_r2)
model_bov_r2_aic <- stepAIC(
  model_bov_r2,
  scope = list(lower = paste0("~`", paste0(magyarazo_rutinmodel, collapse = "`+`"), "`"))
)

summary(model_bov_r2_aic)
confint(model_bov_r2_aic)
AIC(model_bov_r2_aic)

valasztott_valtozok2 <- attr(terms(model_bov_r2_aic), "term.labels")

# Results table (AIC_round2.csv — saving commented out, as originally)
coefs <- summary(model_bov_r2_aic)$coefficients
confint_vals <- confint(model_bov_r2_aic)
results_df <- data.frame(
  Variable = rownames(coefs),
  Estimate = round(coefs[, "Estimate"], 3),
  Std.Error = round(coefs[, "Std. Error"], 3),
  z.value = round(coefs[, "z value"], 3),
  p.value = sapply(coefs[, "Pr(>|z|)"], format_p),
  CI.Lower = round(confint_vals[, 1], 3),
  CI.Upper = round(confint_vals[, 2], 3),
  stringsAsFactors = FALSE
)
# write.csv2(results_df, "AIC_round2.csv", row.names = FALSE)

# Routine model on the r2 patient set (id-based matching)
adatok3 <- adatok_rutin[adatok_rutin[[id_variable]] %in% adatok_bovitett_r2[[id_variable]], ]
model_rutin_sub2 <- glm(build_formula(target_variable, magyarazo_rutinmodel),
                        family = "binomial", data = adatok3)
model_rutin_sub2_aic <- stepAIC(model_rutin_sub2)

summary(model_rutin_sub2)
confint(model_rutin_sub2)
AIC(model_rutin_sub2)
summary(model_rutin_sub2_aic)

magyarazo_bovitett_r2

try(plot.roc(adatok_bovitett_r2$CP_binary, predict(model_bov_r2_aic), print.auc = TRUE), silent = TRUE)


# 5. Interactions + ROC.jpg ----------------------------------------------------
# Requires: Setup + Helpers + section 4 (adatok_bovitett_r2).

model_bov_int <- glm(
  build_formula(target_variable, magyarazo_bovitett_r2, interactions = interakcios_tagok_int),
  family = "binomial", data = adatok_bovitett_r2
)
model_bov_int_aic <- stepAIC(
  model_bov_int,
  scope = list(lower = paste0("~`", paste0(magyarazo_rutinmodel, collapse = "`+`"),
                              "`+`Ethnicity`*`Maternal BMI at OGTT", "`"))
)

summary(model_bov_int_aic)
confint(model_bov_int_aic)
AIC(model_bov_int_aic)

try(plot.roc(adatok_bovitett_r2$CP_binary, predict(model_bov_int_aic), print.auc = TRUE), silent = TRUE)

roc_obj <- roc(adatok_bovitett_r2$CP_binary, predict(model_bov_int_aic))
opt <- coords(roc_obj, "best", ret = c("threshold", "sensitivity", "specificity"))

# Saving ROC.jpg with robust device handling (on.exit dev.off)
save_base_plot("ROC.jpg", {
  plot(
    roc_obj,
    col = "#377eb8",
    lwd = 3,
    print.auc = TRUE,
    print.auc.cex = 1.2,
    legacy.axes = TRUE,
    main = "ROC curve and AUC value",
    grid = TRUE
  )
  abline(a = 0, b = 1, lty = 2, col = "gray")
})

vegso_valtozok <- c(valasztott_valtozok,
                    "`Cesarean Section`*`1st time pregnant`",
                    "`Cesarean Section`:b_CordPGC_mmol.L")


# 6. Bayesian model — load + ROC -----------------------------------------------
# Requires: section 4 (adatok_bovitett_r2) + imputalt.rData in the working directory.
# The model was produced by JointAI (see the section 9 "JointAI generation").

if (!file.exists("imputalt.rData")) {
  stop("Missing imputalt.rData — required for the Bayesian section (see JointAI generation).")
}
load("imputalt.rData")

summary(modelimput)
sum_table <- summary(modelimput, type = "stats", p.values = TRUE)

traceplot(modelimput, col = c('#d4af37', '#460E1B', '#D10E3B'), subset = c(1:8), ncol = 4)

# --- Bayes ROC (DETECTED FAULTY BLOCK, with a safeguard) ---
# The original manual design matrix (fixed [,1:14] indices + manual dummies) fails in the
# original TOO: the name in the deletion list ("Neonatalsexfemale") does not match the actual
# factor column name, so a factor stays in and `X %*% beta_hat` blows up.
# Since this block has no saved output (only an on-screen ROC), tryCatch guards it,
# so the pipeline does not stop. (steps.txt: H6 — to be fixed separately if the figure is needed.)
tryCatch({
model_vars <- all.vars(formula(modelimput))
predictor_vars <- model_vars[model_vars != as.character(formula(modelimput)[[2]])]
adatok_pred_bayes <- adatok_bovitett_r2[, 1:14]
colnames(adatok_pred_bayes)[1:14] <- predictor_vars

adatok_pred_bayes$EthnicityBlack <- as.numeric(adatok_pred_bayes$Ethnicity == "Black")
adatok_pred_bayes$EthnicityAsian <- as.numeric(adatok_pred_bayes$Ethnicity == "Asian")
adatok_pred_bayes$EthnicityHispanic <- as.numeric(adatok_pred_bayes$Ethnicity == "Hispanic")
adatok_pred_bayes$firsttimepregnantTRUE <- as.numeric(adatok_pred_bayes$firsttimepregnant == TRUE)
adatok_pred_bayes$CesareanSection1 <- as.numeric(adatok_pred_bayes$CesareanSection == 1)
adatok_pred_bayes$Neonatalsexfemale2 <- as.numeric(adatok_pred_bayes$Neonatalsex %in% c("female", 2))

adatok_pred_bayes$`CesareanSection1:bCordPGCmg` <-
  adatok_pred_bayes$CesareanSection1 * adatok_pred_bayes$bCordPGCmg
adatok_pred_bayes$`firsttimepregnantTRUE:CesareanSection1` <-
  adatok_pred_bayes$firsttimepregnantTRUE * adatok_pred_bayes$CesareanSection1
adatok_pred_bayes$`MaternalBMIatOGTT:EthnicityBlack` <-
  adatok_pred_bayes$MaternalBMIatOGTT * adatok_pred_bayes$EthnicityBlack
adatok_pred_bayes$`MaternalBMIatOGTT:EthnicityAsian` <-
  adatok_pred_bayes$MaternalBMIatOGTT * adatok_pred_bayes$EthnicityAsian
adatok_pred_bayes$`MaternalBMIatOGTT:EthnicityHispanic` <-
  adatok_pred_bayes$MaternalBMIatOGTT * adatok_pred_bayes$EthnicityHispanic
adatok_pred_bayes$MaternalGDM2.17439608395602 <-
  as.numeric(as.numeric(adatok_pred_bayes$MaternalGDM) > 1)

adatok_pred_bayes <- adatok_pred_bayes[, !(colnames(adatok_pred_bayes) %in% c(
  "Ethnicity", "firsttimepregnant", "CesareanSection", "Neonatalsexfemale", "MaternalGDM"
))]

post_samples <- summary(modelimput, pars = "beta")$res$CPbinary$regcoef
beta_hat <- post_samples[, "Mean"]
X <- as.matrix(adatok_pred_bayes)
X <- cbind(Intercept = 1, X)
linpred <- X %*% beta_hat
pred_vals_bay <- plogis(linpred)

roc_obj_bay <- roc(adatok_bovitett_r2$CP_binary, pred_vals_bay)

plot(
  roc_obj_bay,
  col = "#377eb8",
  lwd = 3,
  print.auc = TRUE,
  print.auc.cex = 1.2,
  legacy.axes = FALSE,
  main = "ROC curve and AUC value",
  grid = TRUE
)
abline(a = 0, b = 1, lty = 2, col = "gray")
opt <- coords(roc_obj_bay, "best", ret = c("threshold", "sensitivity", "specificity"))
points(opt["specificity"], opt["sensitivity"], pch = 19, col = "red", cex = 1.5)
text(
  x = opt["specificity"], y = opt["sensitivity"],
  labels = paste0("Thresh=", round(opt["threshold"], 3)), pos = 4, col = "red"
)
}, error = function(e)
  message("[Bayes ROC skipped] the manual design matrix is faulty (original error): ", conditionMessage(e)))


# 7. Forestplot — GLM vs Bayes -------------------------------------------------
# Requires: section 5 (model_bov_int_aic) + section 6 (modelimput).
# DETECTED FRAGILE BLOCK: the number of coefficients of the saved Bayesian model
# (imputalt.rData) may differ from the fresh stepAIC model (22 vs 23 on the current data), so
# `bayes_coefs$term <- glm_coefs$term` blows up. tryCatch guards the pipeline. If the
# Bayes vs GLM comparison figure is needed, imputalt.rData must be aligned to the fresh model.
tryCatch({
glm_coefs <- broom::tidy(model_bov_int_aic, conf.int = TRUE)
post_samples <- summary(modelimput, pars = "beta")$res$CPbinary$regcoef

bayes_coefs <- data.frame(
  term = rownames(post_samples),
  estimate = post_samples[, "Mean"],
  conf.low = post_samples[, "2.5%"],
  conf.high = post_samples[, "97.5%"],
  tail_prob = post_samples[, "tail-prob."]
)

bayes_coefs$term <- glm_coefs$term
bayes_coefs$tail_prob <- sprintf("%.2e", bayes_coefs$tail_prob)
bayes_coefs$tail_prob[bayes_coefs$tail_prob == "0.00e+00"] <- "<1e-03"
glm_coefs$p.value <- sprintf("%.2e", glm_coefs$p.value)

glm_coefs <- glm_coefs %>%
  select(term, estimate, conf.low = conf.low, conf.high = conf.high, p.value) %>%
  mutate(Model = "Complete cases model")

bayes_coefs <- bayes_coefs %>%
  rename(p.value = tail_prob) %>%
  mutate(Model = "Bayesian imputed")

coefs_combined <- bind_rows(glm_coefs, bayes_coefs)

glm_plot <- glm_coefs %>%
  select(term, est_glm = estimate, low_glm = conf.low, upp_glm = conf.high, p_glm = p.value)
bayes_plot <- bayes_coefs %>%
  select(term, est_bayes = estimate, low_bayes = conf.low, upp_bayes = conf.high, p_bayes = p.value)

plot_data <- left_join(glm_plot, bayes_plot, by = "term") %>%
  mutate(variable = term, var_label = variable, model1 = "GLM", model2 = "Bayes")

plot_data$`Complete case estimates` <- paste0(
  round(exp(plot_data$est_glm), 2), " (", round(exp(plot_data$low_glm), 2),
  " - ", round(exp(plot_data$upp_glm), 2), ")"
)
plot_data$`Bayesian mean estimates` <- paste0(
  round(exp(plot_data$est_bayes), 2), " (", round(exp(plot_data$low_bayes), 2),
  " - ", round(exp(plot_data$upp_bayes), 2), ")"
)
plot_data$`OR (CIs* 2.5%-97.5%)` <- paste(rep(" ", 40), collapse = " ")

colnames(plot_data)[colnames(plot_data) == "p_glm"] <- "Complete case p"
colnames(plot_data)[colnames(plot_data) == "p_bayes"] <- "Imputed p"
colnames(plot_data)[colnames(plot_data) == "variable"] <- "Variables"
plot_data$Variables <- c(
  "(Intercept)", "Maternal age at OGTT", "Maternal BMI at OGTT", "Ethnicity - Afro-Caribbean",
  "Ethnicity - Thai", "Ethnicity - Hispanic", "1st time delivery", "Cesarean Delivery",
  "Maternal GDM", "Age of gestation at delivery", "Birthweight", "Neonatal sex (female)",
  "Maternal HbA1c at OGTT", "Neonatal head circumference", "Neonatal length", "HOMA2-IR at OGTT",
  "Maternal weight gain (till OGTT)", "Neonatal Cord PG concentration",
  "Cesarean Delivery * Neonatal Cord PG concentration", "1st time delivery * Cesarean Delivery",
  "Maternal BMI at OGTT * Ethnicity - Afro-Caribbean", "Maternal BMI at OGTT * Ethnicity - Thai",
  "Maternal BMI at OGTT * Ethnicity - Hispanic"
)

p <- forest(
  data = plot_data[, c(10, 14, 15, 16, 5, 9)],
  est = list(exp(plot_data$est_glm), exp(plot_data$est_bayes)),
  lower = list(exp(plot_data$low_glm), exp(plot_data$low_bayes)),
  upper = list(exp(plot_data$upp_glm), exp(plot_data$upp_bayes)),
  ci_column = 4,
  ref_line = 1,
  arrow_lab = c("Lower risk", "Higher risk"),
  xlim = c(0, 5),
  ticks_at = c(0.5, 1, 2, 3),
  theme = tm_double
)
plot(p)
save_forest(p, "forestplot_models.jpg", width = 5600, height = 4000)
}, error = function(e) message("[forestplot_models skipped]: ", conditionMessage(e)))


# 8. Birth type (joint variable) + combined Bayesian forestplot -------------------------

# --- GLM (complete cases) ---
model_csaszar <- birthtype_model(adatok_master, target_variable)
summary(model_csaszar)

# (labels_birthtype / labels_birthtype_bayes are defined in the Helpers block)
co <- broom::tidy(model_csaszar, conf.int = TRUE)
or_forest_single(
  co$estimate, co$conf.low, co$conf.high, sprintf("%.2e", co$p.value),
  labels_birthtype, "forestplot_csaszar_int.jpg", fallback_labels = co$term
)

# --- ROC for the birth-type GLM (the GLM side of supp_fig7_2model) -> supp_fig7_roc.jpg ---
roc_csaszar <- roc(model_csaszar$y, predict(model_csaszar))
save_base_plot("supp_fig7_roc.jpg", {
  plot(roc_csaszar, col = "#377eb8", lwd = 3, print.auc = TRUE, print.auc.cex = 1.2,
       legacy.axes = TRUE, main = "ROC curve - birth type GLM model", grid = TRUE)
  abline(a = 0, b = 1, lty = 2, col = "gray")
})

# --- Bayes (combined interaction variable) — DETECTED FRAGILE BLOCK ---
# (the saved Bayesian model and the fresh label/term structure may differ; tryCatch guards it)
tryCatch({
if (!file.exists("imputalt_osszevont.rData"))
  stop("Missing imputalt_osszevont.rData — required for the combined Bayesian section.")
load("imputalt_osszevont.rData")

summary(modelimput)
sum_table <- summary(modelimput, type = "stats", p.values = TRUE)

post_samples <- summary(modelimput, pars = "beta")$res$CPbinary$regcoef
bayes_coefs <- data.frame(
  term = rownames(post_samples),
  estimate = post_samples[, "Mean"],
  conf.low = post_samples[, "2.5%"],
  conf.high = post_samples[, "97.5%"],
  tail_prob = post_samples[, "tail-prob."]
)
bayes_coefs$tail_prob <- sprintf("%.2e", bayes_coefs$tail_prob)
bayes_coefs$tail_prob[bayes_coefs$tail_prob == "0.00e+00"] <- "<1e-03"

or_forest_single(
  bayes_coefs$estimate, bayes_coefs$conf.low, bayes_coefs$conf.high, bayes_coefs$tail_prob,
  labels_birthtype, "forestplot_bayes_ov_int_models.jpg",
  est_header = "Bayesian mean estimates", p_header = "Tail-probability",
  xlim = c(0, 6), ticks = c(0.5, 1, 2, 3, 4, 5), width = 4200, height = 3000,
  fallback_labels = bayes_coefs$term
)
}, error = function(e) message("[forestplot_bayes_ov skipped]: ", conditionMessage(e)))


# 9. JointAI generation — (RE)GENERATION of the Bayesian models ----------------
# Requires: Setup + Helpers. Runs ONLY if REGENERATE_BAYES == TRUE (settable at the
# start of Setup). glm_imp (JointAI) is SLOW (several minutes), and it OVERWRITES the
# imputalt.rData / imputalt_osszevont.rData files. It does not run by default; sections 6 and 8
# read in the ready-made .rData.
# NOTE: JointAI requires syntactic names (without spaces/backticks), so we transform the
# column and formula names with gsub (as in the original).

if (REGENERATE_BAYES) {

  ## 9a. Main interaction model -> imputalt.rData
  magyarazo_bovitett_int <- magyarazo_bovitett_r2
  interakcios_tagok <- c(
    "`Cesarean Section`*`b_CordPGC_mmol.L`",
    "`Cesarean Section`*`1st time pregnant`",
    "`Ethnicity`*`Maternal BMI at OGTT`"
  )
  formula_gen <- paste0(target_variable, "~`", paste0(magyarazo_bovitett_int, collapse = "`+`"),
                        "`+", paste0(interakcios_tagok, collapse = "+"))
  formula2 <- gsub(" ", "", formula_gen)
  formula2 <- gsub("_", "", formula2)
  formula2 <- gsub("1st", "first", formula2)
  formula2 <- gsub("(female)", "female", formula2)
  formula2 <- gsub("\\(|\\)", "", formula2)

  adatok_int <- adatok_master
  adatok_int <- rename_vars(adatok_int, lab_map)
  adatok_int <- adatok_int[, c(magyarazo_bovitett_int, "CP_binary")]
  adatok_int <- standardize_predictors(adatok_int, intersect(magyarazo_bovitett_int, colnames(adatok_int)))
  colnames(adatok_int) <- gsub(" ", "", colnames(adatok_int))
  colnames(adatok_int) <- gsub("_", "", colnames(adatok_int))
  colnames(adatok_int) <- gsub("1st", "first", colnames(adatok_int))
  colnames(adatok_int) <- gsub("\\(|\\)", "", colnames(adatok_int))

  modelimput <- glm_imp(as.formula(formula2), family = "binomial", data = adatok_int,
                        n.adapt = 300, n.iter = 5300, progress.bar = 'none', thin = 5,
                        seed = 2020, verbose = TRUE)
  save(modelimput, file = "imputalt.rData")

  ## 9b. Combined "Birth type" model -> imputalt_osszevont.rData
  # The "Birth type" (Cesarean × 1st-time -> 4 categories) + Cesarean:CordPG interaction,
  # with the logic of birthtype_model(), then converted to JointAI-compatible names.
  ov_vars <- magyarazo_bovitett_r2
  adatok_ov <- adatok_master[, c(ov_vars, "CP_binary")]
  adatok_ov <- na.omit(adatok_ov)
  adatok_ov <- standardize_predictors(adatok_ov, ov_vars)
  adatok_ov$`Birth type` <- factor(
    dplyr::case_when(
      adatok_ov$`Cesarean Section` == 0 & adatok_ov$`1st time pregnant` == TRUE  ~ "First born, vaginal delivery",
      adatok_ov$`Cesarean Section` == 1 & adatok_ov$`1st time pregnant` == TRUE  ~ "First born, Cesarean section",
      adatok_ov$`Cesarean Section` == 0 & adatok_ov$`1st time pregnant` == FALSE ~ "Non-first born, vaginal delivery",
      TRUE ~ "Non-first born, Cesarean section"
    ),
    levels = c("First born, vaginal delivery", "Non-first born, vaginal delivery",
               "First born, Cesarean section", "Non-first born, Cesarean section")
  )

  ov_terms <- c(setdiff(ov_vars, c("Cesarean Section", "1st time pregnant")),
                "Birth type", "Cesarean Section`:`b_CordPGC_mmol.L")
  formula_ov <- paste0("CP_binary", "~`", paste0(ov_terms, collapse = "`+`"), "`")
  formula_ov <- gsub(" ", "", formula_ov)
  formula_ov <- gsub("_", "", formula_ov)
  formula_ov <- gsub("1st", "first", formula_ov)
  formula_ov <- gsub("(female)", "female", formula_ov)
  formula_ov <- gsub("\\(|\\)", "", formula_ov)

  colnames(adatok_ov) <- gsub(" ", "", colnames(adatok_ov))
  colnames(adatok_ov) <- gsub("_", "", colnames(adatok_ov))
  colnames(adatok_ov) <- gsub("1st", "first", colnames(adatok_ov))
  colnames(adatok_ov) <- gsub("\\(|\\)", "", colnames(adatok_ov))

  modelimput <- glm_imp(as.formula(formula_ov), family = "binomial", data = adatok_ov,
                        n.adapt = 300, n.iter = 5300, progress.bar = 'none', thin = 5,
                        seed = 2020, verbose = TRUE)
  save(modelimput, file = "imputalt_osszevont.rData")

  rm(modelimput)  # do not pollute the downstream sections 6/8 (they reload their own)
}


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


# 12. Additional analyses for review -------------------------------------------

# 12a. Thai (Asian) population distribution by known birth order
adatok_asian_thai <- adatok[adatok$Ethnicity == "Asian", ]
adatok_asian_thai$known_birth_order <-
  ifelse(is.na(adatok_asian_thai$`1st time pregnant`), 0, 1)
adatok_asian_thai$known_birth_order <- as.factor(adatok_asian_thai$known_birth_order)

vars <- setdiff(names(adatok_asian_thai), "known_birth_order")
table1 <- CreateTableOne(
  vars = vars,
  strata = "known_birth_order",
  data = adatok_asian_thai,
  factorVars = vars[sapply(adatok_asian_thai[vars], is.factor)]
)
print(table1, showAllLevels = TRUE, quote = FALSE, noSpaces = TRUE)

table1_mat <- print(table1, showAllLevels = TRUE, quote = FALSE,
                    noSpaces = TRUE, printToggle = FALSE)
table1_df <- as.data.frame(table1_mat)
table1_df <- cbind(Variable = rownames(table1_df), table1_df)
rownames(table1_df) <- NULL

ft <- flextable(table1_df)
ft <- autofit(ft)
doc <- read_docx()
doc <- body_add_par(doc, "Table 1 – Descriptive statistics", style = "heading 1")
doc <- body_add_flextable(doc, ft)
print(doc, target = "table1_known_birth_order.docx")

model_asian_bo <- glm(
  known_birth_order ~ `Age of gestation at delivery` +
    `Maternal age at OGTT` + `Maternal BMI at OGTT` + Birthweight +
    m_HbA1c_percent + `Cesarean Section` + `m_Smoker.Yes.1.No.0.`,
  data = adatok_asian_thai, family = binomial
)
summary(model_asian_bo)

# 12b. Missingness model (birth order missing ~ CP + covariates)
adatok_asian_thai$birth_order_missing <- as.integer(adatok_asian_thai$known_birth_order == 0)
adatok_asian_thai$Birthweight_100g <- adatok_asian_thai$Birthweight / 100

model_missing <- glm(
  birth_order_missing ~ CP_binary + `Maternal age at OGTT` +
    `Maternal BMI at OGTT` + Birthweight_100g,
  family = binomial, data = adatok_asian_thai
)
summary(model_missing)

tab_missing <- tidy(model_missing, conf.int = TRUE, exponentiate = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    Variable = c(
      "CP positive", "Maternal age at OGTT (per year)",
      "Maternal BMI at OGTT (per kg/m²)", "Birthweight (per 100 g)"
    ),
    OR_CI = sprintf("%.2f (%.2f–%.2f)", estimate, conf.low, conf.high),
    p = signif(p.value, 3)
  ) %>%
  select(Variable, OR_CI, p)
tab_missing

# 12c. Post-hoc power for ancestry
# REQUIRES SECTION 11 (sp_* / model_bov_* ancestry models + adatok_*).
# NOTE: the Asian dataset here — matching the ORIGINAL code — is the data frame of the
# Thai section (12a) above (not the na.omit model data), so that the power value
# stays unchanged. (steps.txt: an intentionally preserved effect of the earlier name collision.)
models <- list(
  EU = model_bov_eu, Black = model_bov_black,
  Hispanic = model_bov_hisp, Asian = model_bov_asian
)
datasets <- list(
  EU = adatok_eu, Black = adatok_black,
  Hispanic = adatok_hisp, Asian = adatok_asian_thai
)

calc_power <- function(model, data) {
  coefs <- summary(model)$coefficients
  beta <- coefs[5, "Estimate"]      # row 5 = exposure
  se <- coefs[5, "Std. Error"]
  OR <- exp(beta)
  n <- nrow(data)
  y_name <- all.vars(formula(model))[1]
  p0 <- mean(data[[y_name]] == 1, na.rm = TRUE)
  p1 <- (OR * p0) / (1 - p0 + OR * p0)
  h <- ES.h(p1, p0)
  power <- pwr.2p.test(h = h, n = n, sig.level = 0.05, alternative = "two.sided")$power
  data.frame(N = n, beta = beta, OR = OR, SE = se,
             outcome_prevalence = p0, p1_est = p1, cohen_h = h, power = power)
}

results <- bind_rows(lapply(names(models), function(g) {
  res <- calc_power(models[[g]], datasets[[g]])
  res$ancestry <- g
  res
}))
results <- results %>%
  select(ancestry, N, beta, OR, SE, outcome_prevalence, p1_est, cohen_h, power)
print(results)


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


# 15. Skinfold models / BMI_prepreg -------------------------------------------
# Requires: Setup + Helpers. "Final variables" = magyarazo_bovitett_r2 + the two main
# interactions (final_interactions, defined in Helpers). The labels are generated from the
# pretty_from_terms() dictionary.

# 15a. Skinfold: the final variables + ALL THREE neonatal skinfolds
#      (flank b_N_FLMn, subscapular b_N_SSMn, triceps b_N_TRMn)
skin_vars <- c(unique(magyarazo_bovitett_r2), "b_N_FLMn", "b_N_SSMn", "b_N_TRMn")
adatok_skin <- na.omit(adatok_master[, c(skin_vars, target_variable)])
adatok_skin <- standardize_predictors(adatok_skin, skin_vars)

model_skinfold <- glm(
  build_formula(target_variable, c(skin_vars, final_interactions)),
  family = "binomial", data = adatok_skin
)
summary(model_skinfold)

co <- broom::tidy(model_skinfold, conf.int = TRUE)
or_forest_single(
  co$estimate, co$conf.low, co$conf.high, sprintf("%.2e", co$p.value),
  pretty_from_terms(co$term), "skinfold_model.jpg", fallback_labels = co$term
)

# 15b. BMI_prepreg: in the final variables, replace "Maternal BMI at OGTT" -> pre_BMI
bmi_vars <- c(setdiff(unique(magyarazo_bovitett_r2), "Maternal BMI at OGTT"), "pre_BMI")
adatok_bmiprepreg <- adatok_master[, c(bmi_vars, target_variable)]
adatok_bmiprepreg$pre_BMI <- as.numeric(as.character(adatok_bmiprepreg$pre_BMI))
adatok_bmiprepreg <- na.omit(adatok_bmiprepreg)
adatok_bmiprepreg <- standardize_predictors(adatok_bmiprepreg, bmi_vars)

model_bmiprepreg <- glm(
  build_formula(target_variable, c(bmi_vars, final_interactions)),
  family = "binomial", data = adatok_bmiprepreg
)
summary(model_bmiprepreg)

co <- broom::tidy(model_bmiprepreg, conf.int = TRUE)
or_forest_single(
  co$estimate, co$conf.low, co$conf.high, sprintf("%.2e", co$p.value),
  pretty_from_terms(co$term), "bmi_prepreg_model.jpg", fallback_labels = co$term
)


# 16. Post-hoc power — total (crude) vs restricted (final model) ---------------
# Requires: Setup (adatok) for TOTAL; section 11 (sp_white/black/asian/hisp) for RESTRICTED.
# Two comparable tables:
#   TOTAL      — the full ancestry sample (where ethnicity + CD + outcome are available), crude CD OR.
#   RESTRICTED — the final-model complete-case set per ancestry, adjusted CD OR.
# The power uses the same pwr.2p.test logic as section 12c.

# Extract the CD (exposure) coefficient BY NAME (robust, not by row index).
cd_coef <- function(model) {
  cf <- summary(model)$coefficients
  rn <- rownames(cf)
  hit <- grep("Cesarean", rn)          # the `Cesarean Section`1 main-effect row
  hit <- hit[!grepl(":", rn[hit])]     # exclude the interaction terms
  cf[hit[1], , drop = FALSE]
}

power_row <- function(model, data, ancestry, y_name = "CP_binary") {
  cc <- cd_coef(model)
  beta <- cc[1, "Estimate"]; se <- cc[1, "Std. Error"]; OR <- exp(beta)
  n <- nrow(data)
  p0 <- mean(data[[y_name]] == 1, na.rm = TRUE)
  p1 <- (OR * p0) / (1 - p0 + OR * p0)
  h <- ES.h(p1, p0)
  power <- pwr.2p.test(h = h, n = n, sig.level = 0.05, alternative = "two.sided")$power
  data.frame(ancestry = ancestry, N = n, beta = beta, OR = OR, SE = se,
             outcome_prevalence = p0, p1_est = p1, cohen_h = h, power = power)
}

anc_levels <- c(EU = "White", Black = "Black", Asian = "Asian", Hispanic = "Hispanic")

# (1) TOTAL: full ancestry sample + crude (unadjusted) CD -> outcome model
power_total <- bind_rows(lapply(names(anc_levels), function(g) {
  d <- adatok[adatok$Ethnicity == anc_levels[[g]], ]
  m <- glm(CP_binary ~ `Cesarean Section`, family = binomial, data = d)
  power_row(m, d, g)
}))

# (2) RESTRICTED: complete-case final ancestry models (section 11 sp_* objects)
sp_list <- list(EU = sp_white, Black = sp_black, Asian = sp_asian, Hispanic = sp_hisp)
power_restricted <- bind_rows(lapply(names(sp_list), function(g) {
  power_row(sp_list[[g]]$model, sp_list[[g]]$data, g)
}))

cat("\n=== Post-hoc power — TOTAL (full ancestry sample, crude CD OR) ===\n")
print(power_total)
cat("\n=== Post-hoc power — RESTRICTED (final-model complete-case, adjusted CD OR) ===\n")
print(power_restricted)

write.csv(power_total, "power_total.csv", row.names = FALSE)
write.csv(power_restricted, "power_restricted.csv", row.names = FALSE)


# 17. Supplementary Fig. 7 — two-model (GLM + Bayes) birth type ----------------
# Requires: section 8 (model_csaszar) + imputalt_osszevont.rData.
# Same STYLE/SIZE/RESOLUTION as forestplot_models.jpg (Fig 1): complete-case GLM vs imputed
# Bayes birth-type model in two columns (tm_double). NEW file (supp_fig7_2model.jpg); the
# forestplot_csaszar_int.jpg is NOT overwritten.
tryCatch({
  e_ov <- new.env(); load("imputalt_osszevont.rData", envir = e_ov)
  modelimput_ov <- e_ov$modelimput

  glm_coefs <- broom::tidy(model_csaszar, conf.int = TRUE)
  post_samples <- summary(modelimput_ov, pars = "beta")$res$CPbinary$regcoef

  bayes_coefs <- data.frame(
    term = rownames(post_samples),
    estimate = post_samples[, "Mean"],
    conf.low = post_samples[, "2.5%"],
    conf.high = post_samples[, "97.5%"],
    tail_prob = post_samples[, "tail-prob."]
  )
  if (nrow(bayes_coefs) != nrow(glm_coefs))
    stop(sprintf("GLM (%d) and Bayes (%d) coefficient counts differ.", nrow(glm_coefs), nrow(bayes_coefs)))

  bayes_coefs$term <- glm_coefs$term
  bayes_coefs$tail_prob <- sprintf("%.2e", bayes_coefs$tail_prob)
  bayes_coefs$tail_prob[bayes_coefs$tail_prob == "0.00e+00"] <- "<1e-03"
  glm_coefs$p.value <- sprintf("%.2e", glm_coefs$p.value)

  glm_plot <- glm_coefs %>%
    select(term, est_glm = estimate, low_glm = conf.low, upp_glm = conf.high, p_glm = p.value)
  bayes_plot <- bayes_coefs %>%
    select(term, est_bayes = estimate, low_bayes = conf.low, upp_bayes = conf.high, p_bayes = tail_prob)

  plot_data <- left_join(glm_plot, bayes_plot, by = "term")

  plot_data$`Complete case estimates` <- paste0(
    round(exp(plot_data$est_glm), 2), " (", round(exp(plot_data$low_glm), 2),
    " - ", round(exp(plot_data$upp_glm), 2), ")"
  )
  plot_data$`Bayesian mean estimates` <- paste0(
    round(exp(plot_data$est_bayes), 2), " (", round(exp(plot_data$low_bayes), 2),
    " - ", round(exp(plot_data$upp_bayes), 2), ")"
  )
  plot_data$`OR (CIs* 2.5%-97.5%)` <- paste(rep(" ", 40), collapse = " ")
  # Local relabel: "Cesarean Section/section" -> "Cesarean Delivery/delivery" (case-preserving),
  # ONLY on this figure. labels_birthtype itself is unchanged (section 8 figures are not affected).
  plot_data$Variables <- gsub("Cesarean section", "Cesarean delivery",
                              gsub("Cesarean Section", "Cesarean Delivery", labels_birthtype))

  disp <- plot_data[, c("Variables", "Complete case estimates", "Bayesian mean estimates",
                        "OR (CIs* 2.5%-97.5%)", "p_glm", "p_bayes")]
  names(disp)[names(disp) == "p_glm"] <- "Complete case p"
  names(disp)[names(disp) == "p_bayes"] <- "Imputed p"

  p <- forest(
    data = disp,
    est = list(exp(plot_data$est_glm), exp(plot_data$est_bayes)),
    lower = list(exp(plot_data$low_glm), exp(plot_data$low_bayes)),
    upper = list(exp(plot_data$upp_glm), exp(plot_data$upp_bayes)),
    ci_column = 4,
    ref_line = 1,
    arrow_lab = c("Lower risk", "Higher risk"),
    xlim = c(0, 10),
    ticks_at = c(0.5, 1, 2, 4, 6, 8),
    theme = tm_double
  )
  plot(p)
  save_forest(p, "supp_fig7_2model.jpg", width = 5600, height = 4000)
}, error = function(e) message("[supp_fig7_2model skipped]: ", conditionMessage(e)))


# 18. Pre-pregnancy BMI vs BMI at OGTT — correlation ---------------------------
# Requires: Setup (adatok). The model covariate is BMI-at-OGTT; pre_BMI is only used for weight gain.
bmi_ogtt <- as.numeric(as.character(adatok$`Maternal BMI at OGTT`))
bmi_pre  <- as.numeric(as.character(adatok$pre_BMI))
ok <- is.finite(bmi_ogtt) & is.finite(bmi_pre)
d_bmi <- data.frame(ogtt = bmi_ogtt[ok], pre = bmi_pre[ok])

set.seed(2020)
sw <- function(x) { xs <- if (length(x) > 5000) sample(x, 5000) else x; shapiro.test(xs) }
sw_ogtt <- sw(d_bmi$ogtt); sw_pre <- sw(d_bmi$pre)
pear  <- cor.test(d_bmi$ogtt, d_bmi$pre, method = "pearson")
spear <- suppressWarnings(cor.test(d_bmi$ogtt, d_bmi$pre, method = "spearman"))

sink("bmi_correlation.txt")
cat("Pre-pregnancy BMI vs BMI at OGTT — correlation\n")
cat("n =", nrow(d_bmi), "complete pairs\n\n")
cat("Normality (Shapiro-Wilk, <=5000 sampled values):\n")
cat(sprintf("  BMI at OGTT : W = %.4f, p = %s\n", sw_ogtt$statistic, format.pval(sw_ogtt$p.value)))
cat(sprintf("  Pre-preg BMI: W = %.4f, p = %s\n\n", sw_pre$statistic, format.pval(sw_pre$p.value)))
cat(sprintf("Pearson  r   = %.3f (95%% CI %.3f to %.3f), p = %s\n",
            pear$estimate, pear$conf.int[1], pear$conf.int[2], format.pval(pear$p.value)))
cat(sprintf("Spearman rho = %.3f, p = %s\n", spear$estimate, format.pval(spear$p.value)))
sink()
cat(readLines("bmi_correlation.txt"), sep = "\n"); cat("\n")
