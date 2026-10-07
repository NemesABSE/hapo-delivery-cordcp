## 2026-10-06: Fig. 2a regenerated from point2_within_stratum_contrasts/point2_fig2a_annotated.R with two display-only label changes:
## the interaction row reads "1st time delivery * Cesarean delivery*,**" (comma raised to the level of the asterisks, as in
## suppfig_stars_1006) and "Cesarean Delivery" is written "Cesarean delivery". Same fit and same Bayesian draws.
## Run from this folder: Rscript fig2a_stars.R > fig2a_stars.log 2>&1
Sys.setenv(REV2_ROOT = normalizePath("..")); source("../_common/setup.R")
pdf(NULL)
REV <- file.path(SHARE, "L/CP/cikkhez/review")
getB <- function(path) {
  e <- new.env(); load(path, envir = e); m <- get(ls(e)[1], envir = e)
  cl <- m$coef_list$CPbinary; M <- do.call(rbind, lapply(m$MCMC, as.matrix))
  B <- M[, cl$coef, drop = FALSE]; colnames(B) <- cl$varname
  cat(basename(path), ": rows", nrow(m$data), " NA cells", sum(is.na(m$data)), " draws", nrow(B), " betas", ncol(B), "\n"); B
}
norm1 <- function(x) sub("^1sttimepregnantTRUE$", "firsttimepregnantTRUE", gsub("[^A-Za-z0-9]", "", x))
normterm <- function(x) { parts <- strsplit(x, ":")[[1]]; paste(sort(norm1(parts), method = "radix"), collapse = ":") }
match_cols <- function(glm_terms, bnames) {
  gn <- sapply(glm_terms, normterm); bn <- sapply(bnames, normterm)
  bn <- sub("^MaternalGDM[0-9.\\-]+$", "MaternalGDM", bn)                       # JointAI turned GDM into a factor level name
  idx <- match(gn, bn)
  if (any(is.na(idx))) stop("unmatched glm terms: ", paste(glm_terms[is.na(idx)], collapse = " | "), "\nJointAI: ", paste(bnames, collapse = " | "))
  stopifnot(!anyDuplicated(idx)); idx
}
tpstr <- function(z) sprintf("%.2e", 2 * pnorm(-abs(mean(z) / sd(z))))
ortxt <- function(e, l, u) paste0(round(exp(e), 2), " (", round(exp(l), 2), " - ", round(exp(u), 2), ")")

mA <- fit_final(); gA <- broom::tidy(mA, conf.int = TRUE)
cat("Final model: n =", nobs(mA), " AUC =", auc_of(mA, d_cc_s), " terms =", nrow(gA), "\n"); print(gA$term)
BA <- getB(file.path(REV, "imputalt.rData")); cat("JointAI varnames:\n"); print(colnames(BA))
idx <- match_cols(gA$term, colnames(BA)); blist <- lapply(idx, function(i) BA[, i])
print(data.frame(glm = gA$term, jointai = colnames(BA)[idx]))

labA <- c("(Intercept)", "Maternal age at OGTT", "Maternal BMI at OGTT", "Ethnicity - Afro-Caribbean",
  "Ethnicity - Thai", "Ethnicity - Hispanic", "1st time delivery*", "Cesarean delivery**", "Maternal GDM",
  "Age of gestation at delivery", "Birthweight", "Neonatal sex (female)", "Maternal HbA1c at OGTT",
  "Neonatal head circumference", "Neonatal length", "HOMA2-IR at OGTT", "Maternal weight gain (till OGTT)",
  "Neonatal Cord PG concentration", "1st time delivery * Cesarean delivery*,**",
  "Cesarean delivery * Neonatal Cord PG concentration", "Maternal BMI at OGTT * Ethnicity - Afro-Caribbean",
  "Maternal BMI at OGTT * Ethnicity - Thai", "Maternal BMI at OGTT * Ethnicity - Hispanic")
## label order check against the glm terms
chk <- c("Intercept", "age", "BMI", "Black", "Asian", "Hispanic", "1st time", "Cesarean", "GDM", "gestation", "Birthweight",
         "sex", "HbA1c", "HCMn", "LNGMn", "HOMA2", "Wtgain|WtMgain", "CordPGC", "1st time.*Cesarean|Cesarean.*1st time",
         "CordPGC.*Cesarean|Cesarean.*CordPGC", "Black", "Asian", "Hispanic")
stopifnot(all(mapply(function(p, t) grepl(p, t), chk, gA$term)))

best <- sapply(blist, mean); blow <- sapply(blist, quantile, .025); bupp <- sapply(blist, quantile, .975)
disp <- data.frame(Variables = labA, `Complete case estimates` = ortxt(gA$estimate, gA$conf.low, gA$conf.high),
  `Bayesian estimates` = ortxt(best, blow, bupp), `OR (CIs* 2.5%-97.5%)` = paste(rep(" ", 40), collapse = " "),
  `Complete case p` = sprintf("%.2e", gA$p.value), `Imputed p` = sapply(blist, tpstr), check.names = FALSE, stringsAsFactors = FALSE)
FOOT <- "* 1st vs non-first delivery among vaginal deliveries.    ** Cesarean vs vaginal delivery among non-firstborns."
p <- forest(data = disp, est = list(exp(gA$estimate), exp(best)), lower = list(exp(gA$conf.low), exp(blow)),
            upper = list(exp(gA$conf.high), exp(bupp)), ci_column = 4, ref_line = 1,
            arrow_lab = c("Lower risk", "Higher risk"), xlim = c(0, 5), ticks_at = c(0.5, 1, 2, 3),
            theme = tm_double)   # no footnote: the * / ** are explained in the figure legend
raise_comma <- function(p) {
  i <- which(sapply(p$grobs, function(g) inherits(g, "text") && grepl("\\*,\\*\\*$", g$label))); stopifnot(length(i) == 1); g <- p$grobs[[i]]
  a <- sub(",\\*\\*$", "", g$label); y <- g$y - 0.5 * grid::stringHeight(g$label)
  tg <- function(lab, dx, dy = 0) grid::textGrob(lab, x = g$x + dx, y = y + grid::unit(dy, "char"), hjust = 0, vjust = 0)
  p$grobs[[i]] <- grid::gTree(children = grid::gList(tg(a, grid::unit(0, "pt")), tg(",", grid::stringWidth(a), 0.38), tg("**", grid::stringWidth(paste0(a, ",")))), gp = g$gp, vp = g$vp); p }
save_forest(raise_comma(p), "Figure_2a_final_model.jpg", width = 5600, height = 4000)
out <- data.frame(term = gA$term, label = labA, glm_OR = exp(gA$estimate), glm_lo = exp(gA$conf.low), glm_hi = exp(gA$conf.high),
                  glm_p = gA$p.value, bayes_OR = exp(best), bayes_lo = exp(blow), bayes_hi = exp(bupp), bayes_tailp = as.numeric(sapply(blist, tpstr)))
write.csv(out, "Figure_2a_coefficients.csv", row.names = FALSE)
cat("\nKey rows (published: CD 2.23 / 2.74; 1st 0.51 / 0.59; 1st*CD 1.96 / 1.77; CD*PG 1.25 / 1.23):\n")
print(disp[c(7, 8, 19, 20), 1:3], row.names = FALSE)
cat("DONE\n")
