## Shared setup for the second-revision analyses.
## Runs the Setup + Helpers blocks (lines 1-510) of R/02_supplementary_forests.r, so data preparation, variable
## lists, standardization and the AIC-selected variable set (R/stepAIC_results.rds) are those of the main analysis.
## Data: HAPO_DATA (see config.R). Saved JointAI fits, where a script needs them: HAPO_BAYES.
## Packages: MASS, broom, dplyr, grid, pROC, forestploter (loaded by the pipeline); base R otherwise.
ORIG_WD <- getwd()   # outputs go next to the calling script
ROOT <- normalizePath(Sys.getenv("REV2_ROOT", unset = file.path(getwd(), "..")))
REPO <- normalizePath(file.path(ROOT, "..", ".."))
source(file.path(REPO, "config.R"))
hapo_bayes_dir <- function() {
  d <- Sys.getenv("HAPO_BAYES", unset = NA)
  if (is.na(d) || !dir.exists(d)) stop("HAPO_BAYES is not set. It must point at the directory holding the saved JointAI fits (see README).", call. = FALSE)
  d
}
src <- readLines(file.path(REPO, "R", "02_supplementary_forests.r"))[1:510]
src <- src[!grepl('^\\.libPaths\\(', src)]
src <- sub('^setwd\\(normalizePath\\("figure-relabel/run"\\)\\)', sprintf('setwd("%s")', file.path(REPO, "R")), src)
eval(parse(text = src), envir = globalenv())
stopifnot(exists("magyarazo_bovitett_r2"), exists("final_interactions"))
magyarazo_bovitett_r2 <- unique(magyarazo_bovitett_r2)   # duplicates collapse in a formula anyway
int3 <- c(final_interactions, "Ethnicity`*`Maternal BMI at OGTT")   # the 3 final-model interactions

orci <- function(m) {
  co <- broom::tidy(m, conf.int = TRUE)
  data.frame(term = co$term, OR = exp(co$estimate), lo = exp(co$conf.low), hi = exp(co$conf.high),
             p = co$p.value, stringsAsFactors = FALSE)
}
fmt <- function(x) sprintf("%.2f (%.2f to %.2f)", x[1], x[2], x[3])
contrast <- function(m, L) {            # delta-method contrast of log-ORs, Wald 95% CI
  b <- coef(m); V <- vcov(m); L <- L[names(b)]; L[is.na(L)] <- 0
  est <- sum(L * b); se <- sqrt(as.numeric(t(L) %*% V %*% L))
  c(exp(est), exp(est - 1.96 * se), exp(est + 1.96 * se), p = 2 * pnorm(-abs(est / se)))
}
find_term <- function(m, pattern) { nm <- names(coef(m)); nm[grepl(pattern, nm)] }
birthtype_contrasts <- function(m) {
  nm <- names(coef(m)); z <- setNames(rep(0, length(nm)), nm)
  b_first <- find_term(m, "^`1st time pregnant`TRUE$"); b_cd <- find_term(m, "^`Cesarean Section`1$")
  b_int <- find_term(m, "Cesarean Section`1:`1st time pregnant`TRUE$|1st time pregnant`TRUE:`Cesarean Section`1$")
  L1 <- z; L1[b_first] <- -1
  L2 <- z; L2[b_cd] <- 1; L2[b_int] <- 1
  L3 <- z; L3[b_cd] <- 1; L3[b_first] <- -1
  L4 <- z; L4[b_cd] <- 1
  rbind(`Non-first born, vaginal (ref first vaginal)` = contrast(m, L1),
        `First born, cesarean (ref first vaginal)`    = contrast(m, L2),
        `Non-first born, cesarean (ref first vaginal)`= contrast(m, L3),
        `CD vs vaginal within non-firstborns`         = contrast(m, L4),
        `CD vs vaginal within firstborns`             = contrast(m, L2))
}
print_contrasts <- function(cm) for (i in seq_len(nrow(cm)))
  cat(sprintf("  %-48s %s  p = %.2e\n", rownames(cm)[i], fmt(cm[i, 1:3]), cm[i, 4]))
auc_of <- function(m, d) round(as.numeric(pROC::auc(pROC::roc(d[[target_variable]], fitted(m), quiet = TRUE))), 4)

## complete-case set of the final model (n = 3,482), standardized as in the pipeline
d_cc <- adatok_master %>% select(all_of(c(id_variable, magyarazo_bovitett_r2, target_variable))) %>% na.omit()
d_cc_s <- standardize_predictors(d_cc, magyarazo_bovitett_r2)
fit_final <- function(vars = magyarazo_bovitett_r2, ints = int3, data = d_cc_s)
  glm(build_formula(target_variable, c(vars, ints)), family = binomial(), data = data)

## Fig 3 (Birth type) data, CORRECTED: rows with missing birth order are dropped first
## (as in R/01_main_pipeline.r)
bt_data <- function(data_master = adatok_master) {
  vars <- setdiff(magyarazo_bovitett_r2, "1st time pregnant")
  d <- data_master %>% filter(!is.na(`1st time pregnant`)) %>%
    mutate(`Birth type` = factor(case_when(
      `Cesarean Section` == 0 & `1st time pregnant` == TRUE  ~ "First born, vaginal delivery",
      `Cesarean Section` == 1 & `1st time pregnant` == TRUE  ~ "First born, Cesarean section",
      `Cesarean Section` == 0 & `1st time pregnant` == FALSE ~ "Non-first born, vaginal delivery",
      TRUE ~ "Non-first born, Cesarean section"),
      levels = c("First born, vaginal delivery", "Non-first born, vaginal delivery",
                 "First born, Cesarean section", "Non-first born, Cesarean section")))
  vars <- c(vars, "Birth type")
  d <- d %>% select(all_of(c(id_variable, vars, target_variable))) %>% na.omit()
  list(d = standardize_predictors(d, vars), vars = vars)
}
fit_birthtype <- function(with_bmi_ancestry = FALSE, bt = bt_data()) {
  terms <- wrap_terms(c(setdiff(bt$vars, "Cesarean Section"), "Cesarean Section:b_CordPGC_mmol.L"))
  f <- paste(target_variable, "~", paste(terms, collapse = " + "))
  if (with_bmi_ancestry) f <- paste(f, "+ Ethnicity:`Maternal BMI at OGTT`")
  glm(as.formula(f), family = binomial(), data = bt$d)
}
setwd(ORIG_WD)
cat("setup OK: complete-case n =", nrow(d_cc), " | outputs written to", getwd(), "\n")
