# Supplementary Fig. 4, rebuilt as a single five-level delivery-mode-by-birth-order
# variable. Mirrors birthtype_model()/Fig. 2b: the combined categorical exposure plus
# the binary cesarean x cord PG interaction. The sixth cell (repeat cesarean in a
# firstborn) cannot exist, so it is simply not a level and nothing is aliased.

extra <- read.csv(hapo_file("HAPO_Clinical_with_binary_cordcp.csv"))
.i <- match(adatok_cp95[[id_variable]], extra$m_mom_geneva_id)
dt <- extra$m_LD_DelType[.i]
adatok_cp95$DelType3 <- dplyr::case_when(dt %in% c(1,2) ~ "Vaginal",
                                         dt == 3 ~ "Primary", dt == 4 ~ "Repeat")

vars <- setdiff(magyarazo_bovitett_r2, "1st time pregnant")   # keeps `Cesarean Section`
d <- adatok_cp95 %>%
  select(all_of(c(id_variable, vars, "1st time pregnant", "DelType3", "CP_binary"))) %>%
  na.omit()
cat("complete cases:", nrow(d), "\n")
imp <- d$DelType3 == "Repeat" & d$`1st time pregnant` == TRUE
cat("contradictory (repeat cesarean in a firstborn), excluded:", sum(imp), "\n")
d <- d[!imp, ]
cat("n used:", nrow(d), "\n")

LV <- c("First born, vaginal delivery", "Non-first born, vaginal delivery",
        "First born, primary cesarean section", "Non-first born, primary cesarean section",
        "Non-first born, repeat cesarean section")
first <- d$`1st time pregnant` == TRUE
d$`Birth type` <- factor(dplyr::case_when(
  d$DelType3=="Vaginal" & first  ~ LV[1],
  d$DelType3=="Vaginal" & !first ~ LV[2],
  d$DelType3=="Primary" & first  ~ LV[3],
  d$DelType3=="Primary" & !first ~ LV[4],
  d$DelType3=="Repeat"                          ~ LV[5]), levels = LV)
print(table(d$`Birth type`))

mvars <- c(setdiff(vars, "Cesarean Section"), "Birth type")
d <- standardize_predictors(d, mvars)
terms <- c(mvars, "Cesarean Section`:`b_CordPGC_mmol.L")
m <- glm(build_formula("CP_binary", terms), family = "binomial", data = d)
co <- broom::tidy(m, conf.int = TRUE)
cat("\naliased coefficients:", sum(is.na(co$estimate)), "\n")
stopifnot(sum(is.na(co$estimate)) == 0)
write.csv(transform(co, OR=exp(estimate), lo=exp(conf.low), hi=exp(conf.high)),
          "supp_fig4_5level_estimates.csv", row.names = FALSE)

labs <- pretty_from_terms(co$term)
fix <- setNames(LV[-1], paste0("`Birth type`", LV[-1]))
fix <- c(fix, "`Cesarean Section`1:b_CordPGC_mmol.L" =
              "Cesarean delivery * Neonatal Cord PG concentration")
h <- co$term %in% names(fix); labs[h] <- unname(fix[co$term[h]])
print(data.frame(term=co$term, label=labs))

or_forest_single(co$estimate, co$conf.low, co$conf.high,
                 sprintf("%.2e", co$p.value), labs,
                 "SuppFig4_primary_vs_repeat_CD_5level.jpg",
                 xlim = c(0, 11), ticks = c(1, 2, 4, 6, 8, 10),
                 fallback_labels = co$term)

# key contrast: repeat versus primary cesarean among parous women
d2 <- d; d2$`Birth type` <- relevel(d2$`Birth type`, ref = LV[4])
c2 <- broom::tidy(glm(build_formula("CP_binary", terms), family="binomial", data=d2),
                  conf.int = TRUE)
k <- c2[c2$term == paste0("`Birth type`", LV[5]), ]
cat(sprintf("\nRepeat vs PRIMARY cesarean, non-first born: OR %.2f (%.2f to %.2f), P = %.2f\n",
            exp(k$estimate), exp(k$conf.low), exp(k$conf.high), k$p.value))
cat("PLOT DONE\n")
