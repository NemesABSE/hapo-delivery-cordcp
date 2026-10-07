## Reviewer 1, point 6, follow-up (2026-09-09): recursive path analysis behind the DAG, with the
## cord glucose <-> cord C-peptide arrow drawn in ONE direction at a time and standardized path
## coefficients on every substantive arrow; fitted separately in the GDM and non-GDM strata.
## Changes to the DAG relative to Supp_Fig_DAG_point6: (i) no dashed reverse arrow, the two
## directions are two separate models; (ii) birthweight -> cesarean delivery added, so gestational
## age and birthweight are separate nodes (with gestational age -> birthweight); (iii) outcome is
## neonatal hyperinsulinemia (binary, as in the final model; changed 2026-09-11 from log cord C-peptide
## on request), so the arrows into the outcome are odds ratios from a logistic equation.
## Estimation: each endogenous node regressed on its parents plus the final-model covariates
## (base R lm; logistic glm for the binary cesarean and hyperinsulinemia nodes). For a recursive model this is the
## maximum-likelihood path analysis, no SEM package needed. All continuous variables are z-scored
## within stratum, so a slope is the SD change of the child per SD of the parent (per cesarean
## vs vaginal for the binary parent). Interaction terms of the final model are not representable
## as single arrows and are omitted. Run from this folder: Rscript point6_paths.R > point6_paths.log 2>&1
Sys.setenv(REV2_ROOT = normalizePath("..")); source("../_common/setup.R")

## --- data: the 3,482 complete cases of the final model; outcome = binary hyperinsulinemia ----------
d <- d_cc
d$CP_binary <- factor(d$CP_binary, levels = c(0, 1))
cat("complete-case rows:", nrow(d), " | hyperinsulinemia missing:", sum(is.na(d$CP_binary)), "\n")
stopifnot(!anyNA(d$CP_binary))
cat("hyperinsulinemia: "); print(table(d$CP_binary))
cat("GDM strata: "); print(table(d$`Maternal GDM`, useNA = "ifany"))
cat("cesarean by stratum:\n"); print(table(GDM = d$`Maternal GDM`, CD = d$`Cesarean Section`))

CD <- "Cesarean Section"; GA <- "Age of gestation at delivery"; BW <- "Birthweight"
PG <- "b_CordPGC_mmol.L"; CP <- "CP_binary"; GDM <- "Maternal GDM"
covs <- setdiff(magyarazo_bovitett_r2, c(CD, GA, BW, PG, GDM))
cat("covariate block (adjusted in every equation, not drawn):", paste(covs, collapse = "; "), "\n")

bt <- function(v) paste0("`", v, "`")
fit_eq_base <- function(child, parents, data, family = "gaussian", covs) {
  f <- as.formula(paste(bt(child), "~", paste(bt(c(parents, covs)), collapse = " + ")))
  if (family == "gaussian") lm(f, data = data) else glm(f, family = binomial(), data = data)
}
## one row per drawn arrow: coefficient of `parent` in the equation of `child`
edge_rows <- function(m, child, parents, direction, stratum, data) {
  co <- broom::tidy(m, conf.int = TRUE)
  out <- lapply(parents, function(p) {
    r <- co[grepl(paste0("^`?", p, "`?(1|TRUE)?$"), co$term) | co$term == p, ]
    stopifnot(nrow(r) == 1)
    data.frame(direction = direction, stratum = stratum, n = nobs(m),
               equation = if (inherits(m, "glm")) "logistic" else "linear",
               child = child, parent = p, estimate = r$estimate, se = r$std.error,
               ci_lo = r$conf.low, ci_hi = r$conf.high, p = r$p.value,
               scale = if (inherits(m, "glm")) (if (is.factor(data[[p]])) "log-odds, level 1 vs 0 of parent" else "log-odds per SD of parent") else
                       if (is.factor(data[[p]])) "SD of child, level 1 vs 0 of parent" else "SD of child per SD of parent",
               stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

run_stratum <- function(g) {
  ## g = 0 / 1 = GDM stratum (GDM dropped from the covariates); g = "all" = pooled, GDM kept as a covariate
  ds <- if (identical(g, "all")) d else d[d[[GDM]] == g, ]
  covs <- if (identical(g, "all")) c(covs, GDM) else covs
  n_hi <- sum(ds[[CP]] == 1)
  ds <- standardize_predictors(ds, c(GA, BW, PG, covs))         # z-scores within stratum; factors untouched
  fit_eq <- function(child, parents, data, family = "gaussian") fit_eq_base(child, parents, data, family, covs = covs)
  lab <- if (identical(g, "all")) "all" else if (g == 1) "GDM" else "non-GDM"
  cat(sprintf("\n==== stratum %s: n = %d, cesarean = %d, hyperinsulinemia = %d | covariates: %s ====\n",
              lab, nrow(ds), sum(ds[[CD]] == 1), n_hi, paste(covs, collapse = "; ")))
  ## arrows that do not depend on the glucose/C-peptide direction
  m_cd <- fit_eq(CD, BW, ds, "binomial"); m_bw <- fit_eq(BW, GA, ds)
  shared <- rbind(edge_rows(m_cd, CD, BW, "both", lab, ds), edge_rows(m_bw, BW, GA, "both", lab, ds))
  ## direction A: cord glucose -> cord C-peptide
  mA_pg <- fit_eq(PG, c(CD, GA, BW), ds); mA_cp <- fit_eq(CP, c(PG, CD, GA, BW), ds, "binomial")
  A <- rbind(edge_rows(mA_pg, PG, c(CD, GA, BW), "PG->CP", lab, ds),
             edge_rows(mA_cp, CP, c(PG, CD, GA, BW), "PG->CP", lab, ds))
  ## direction B: cord C-peptide -> cord glucose
  mB_cp <- fit_eq(CP, c(CD, GA, BW), ds, "binomial"); mB_pg <- fit_eq(PG, c(CP, CD, GA, BW), ds)
  B <- rbind(edge_rows(mB_cp, CP, c(CD, GA, BW), "CP->PG", lab, ds),
             edge_rows(mB_pg, PG, c(CP, CD, GA, BW), "CP->PG", lab, ds))
  res <- rbind(shared, A, B)
  print(res[, c("direction", "child", "parent", "estimate", "ci_lo", "ci_hi", "p")], digits = 3, row.names = FALSE)

  ## with a logistic outcome the two directions are no longer likelihood-equivalent factorizations and the
  ## direct + indirect = total identity does not hold; report the glucose-adjusted (direct) and the
  ## glucose-omitted (total) cesarean odds ratio instead
  g_ <- function(tab, ch, pa) tab$estimate[tab$child == ch & tab$parent == pa]
  cat(sprintf("CD -> hyperinsulinemia: OR glucose-adjusted (PG->CP model) %.2f | OR glucose omitted (CP->PG model) %.2f\n",
              exp(g_(A, CP, CD)), exp(g_(B, CP, CD))))
  res
}
paths <- rbind(run_stratum(0), run_stratum(1), run_stratum("all"))
write.csv(paths, "point6_path_coefficients.csv", row.names = FALSE)
cat("\nwritten point6_path_coefficients.csv (", nrow(paths), "rows )\nDONE\n")
