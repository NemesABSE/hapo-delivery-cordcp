## 2026-10-06: supplementary forest plots that contain the cesarean delivery x first-time delivery interaction, regenerated
## with the Fig. 2a marks on the row labels: "1st time delivery*", "Cesarean delivery**" and, on the interaction row,
## "1st time delivery * Cesarean delivery*,**" (both marks at the end of the label). Labels are display-only; models are those of pipeline sections 10a
## (Supp Fig 1), 11 (Supp Figs 6 to 9) and 13c (Supp Fig 11c). Also writes, per figure, the cesarean versus vaginal
## delivery contrast within firstborns (main effect + interaction, delta method) for the figure legends.
## Rscript suppfig_stars.R > suppfig_stars.log 2>&1
Sys.setenv(REV2_ROOT = normalizePath("..")); source("../_common/setup.R"); pdf(NULL)
star <- function(l) { l[l == "1st time delivery"] <- "1st time delivery*"; l[l == "Cesarean delivery"] <- "Cesarean delivery**"
  l[l == "1st time delivery * Cesarean delivery"] <- "1st time delivery * Cesarean delivery*,**"; stopifnot(sum(grepl("\\*\\*$", l)) == 2, sum(l == "1st time delivery*") == 1); l }
## The comma of "*,**" is redrawn raised to the level of the asterisks (the label is plain text, so it has no superscript of its own).
raise_comma <- function(p) {
  i <- which(sapply(p$grobs, function(g) inherits(g, "text") && grepl("\\*,\\*\\*$", g$label))); stopifnot(length(i) == 1); g <- p$grobs[[i]]
  a <- sub(",\\*\\*$", "", g$label); y <- g$y - 0.5 * grid::stringHeight(g$label)
  tg <- function(lab, dx, dy = 0) grid::textGrob(lab, x = g$x + dx, y = y + grid::unit(dy, "char"), hjust = 0, vjust = 0)
  p$grobs[[i]] <- grid::gTree(children = grid::gList(tg(a, grid::unit(0, "pt")), tg(",", grid::stringWidth(a), COMMA_RAISE), tg("**", grid::stringWidth(paste0(a, ",")))), gp = g$gp, vp = g$vp); p }
COMMA_RAISE <- 0.38
or_forest_single <- local({ f <- or_forest_single; function(estimate, conf.low, conf.high, p_formatted, labels, file, ..., width = 3600, height = 3000) {
  p <- f(estimate, conf.low, conf.high, p_formatted, labels, file, ..., width = width, height = height); save_forest(raise_comma(p), file, width = width, height = height) } })
res <- list(); coefs <- list()
keep <- function(fig, m, co, labels) {
  b <- coef(m); V <- vcov(m); cd <- "`Cesarean Section`1"
  bi <- find_term(m, "Cesarean Section`1:`1st time pregnant`TRUE$|1st time pregnant`TRUE:`Cesarean Section`1$"); stopifnot(length(bi) == 1)
  est <- b[cd] + b[bi]; se <- sqrt(V[cd, cd] + V[bi, bi] + 2 * V[cd, bi])
  res[[fig]] <<- data.frame(figure = fig, n = nobs(m), CD_nonfirst = exp(b[cd]), CD_first = exp(est), lo = exp(est - 1.96 * se), hi = exp(est + 1.96 * se),
                            p = 2 * pnorm(-abs(est / se)), interaction = exp(b[bi]), row.names = NULL)
  coefs[[fig]] <<- data.frame(figure = fig, term = co$term, label = labels, est = exp(co$estimate), lo = exp(co$conf.low), hi = exp(co$conf.high), p = co$p.value)
}
## Supp Fig 1 (section 10a)
sens_vars <- c(setdiff(magyarazo_bovitett_r2, c("HOMA2_IR", "Maternal GDM")), "m_FCP_nmol_l", "m_FPG_CLC_mmol_l", "m_OneHrPG_CLC_mmol_l", "m_TwoHrPG_CLC_mmol_l")
uj <- adatok %>% select(all_of(id_variable), m_FCP_nmol_l, m_FPG_CLC_mmol_l, m_OneHrPG_CLC_mmol_l, m_TwoHrPG_CLC_mmol_l)
ds <- rename_vars(adatok_master %>% left_join(uj, by = id_variable), lab_map); sens_vars <- rename_vec(sens_vars, lab_map)
ds <- standardize_predictors(na.omit(ds[, c(sens_vars, target_variable)]), sens_vars)
m <- glm(build_formula(target_variable, c(sens_vars, "Cesarean Section`*`1st time pregnant", "Cesarean Section`*`b_CordPGC_mmol.L")), family = "binomial", data = ds)
co <- broom::tidy(m, conf.int = TRUE); L <- star(labels_sens1); stopifnot(length(L) == nrow(co)); keep("SuppFig1", m, co, L)
or_forest_single(co$estimate, co$conf.low, co$conf.high, sprintf("%.2e", co$p.value), L, "SuppFig1_alt_glycemic_covariates_stars.jpg", fallback_labels = co$term)
## Supp Figs 6 to 9 (section 11, as figures/code/regen_ancestry.R)
for (x in list(c("White", "SuppFig6_European"), c("Black", "SuppFig7_AfroCaribbean"), c("Hispanic", "SuppFig8_Hispanic"), c("Asian", "SuppFig9_Thai"))) {
  d <- na.omit(adatok[adatok$Ethnicity == x[1], c(magyarazo_bovitett_r2, target_variable)]); d <- standardize_predictors(d, magyarazo_bovitett_r2)
  v <- c(setdiff(magyarazo_bovitett_r2, "Ethnicity"), "Cesarean Section`*`1st time pregnant", "Cesarean Section`*`b_CordPGC_mmol.L")
  m <- glm(build_formula(target_variable, v), family = "binomial", data = d)
  co <- broom::tidy(m, conf.int = TRUE); L <- star(labels_subpop); stopifnot(length(L) == nrow(co)); keep(x[2], m, co, L)
  or_forest_single(co$estimate, co$conf.low, co$conf.high, sprintf("%.2e", co$p.value), L, paste0(x[2], "_stars.jpg"), fallback_labels = co$term)
}
## Supp Fig 11c (section 13c)
cp_cont <- adatok %>% select(all_of(id_variable), b_CordCP_ug)
dc <- adatok_master %>% left_join(cp_cont, by = id_variable); dc <- na.omit(dc[, c(magyarazo_bovitett_r2, "b_CordCP_ug")])
dc <- dc[dc$b_CordCP_ug > 0, ]; dc <- standardize_predictors(dc, magyarazo_bovitett_r2); dc$logCP <- log(dc$b_CordCP_ug)
m <- lm(build_formula("logCP", c(magyarazo_bovitett_r2, final_interactions)), data = dc)
co <- broom::tidy(m, conf.int = TRUE); L <- star(pretty_from_terms(co$term)); keep("SuppFig11c", m, co, L)
or_forest_single(co$estimate, co$conf.low, co$conf.high, sprintf("%.2e", co$p.value), L, "SuppFig11c_continuous_logCP_stars.jpg",
  est_header = "Fold-change per SD (exp beta)", or_header = "Ratio (95% CI)", xlim = c(0.5, 1.5), ticks = c(0.5, 0.75, 1, 1.25, 1.5),
  width = 5600, height = 3400, arrow_lab = NULL, fallback_labels = co$term)
out <- do.call(rbind, res); print(out, digits = 4, row.names = FALSE)
write.csv(out, "firstborn_contrasts.csv", row.names = FALSE); write.csv(do.call(rbind, coefs), "suppfig_stars_coefficients.csv", row.names = FALSE); cat("DONE\n")
