# Mode of delivery, birth order and neonatal cord blood hyperinsulinemia — analysis code

R code for the analyses reported in *Association of Mode of Delivery and Birth Order with
Hyperinsulinemia in Neonatal Cord Blood* (HAPO substudy).

## Data

The analysis uses the Hyperglycemia and Adverse Pregnancy Outcome (HAPO) substudy dataset,
which is **controlled-access and is not distributed here**. NIH eRA Commons account holders
may request access to the HAPO data through the NIH database of Genotypes and Phenotypes
(dbGaP), study accession `phs000096.v4.p1`, as described at
<https://www.ncbi.nlm.nih.gov/books/NBK570242/>.

The scripts expect three CSV files in one directory:

| File | Contents |
|---|---|
| `Hapo_Clinical_FTO_genotype - extended.csv` | maternal and neonatal clinical variables |
| `HAPO_Clinical_with_binary_cordcp.csv` | cord C-peptide, binary hyperinsulinemia flags, delivery type |
| `Hapo_Clinical_FTO_genotype - extended_with_origPrDel20.csv` | original parity variable, 95th-percentile flag |

Point `HAPO_DATA` at that directory. No path is hard-coded; `config.R` resolves every data
file and fails with an explicit message if the directory is not set.

```r
Sys.setenv(HAPO_DATA = "/path/to/hapo/csv")
source("config.R")
```

## Running the analyses

Run from the repository root, in this order. Scripts 04 and 05 reuse objects created by 02.

| Script | Produces |
|---|---|
| `R/01_main_pipeline.r` | full pipeline: iterative AIC variable selection, final complete-case model, Bayesian joint model, all sensitivity analyses |
| `R/02_supplementary_forests.r` | setup plus the supplementary forest plots (sensitivity, ancestry strata, 95th-percentile outcome, parity) |
| `R/04_suppfig4_five_level.r` | Supplementary Fig. 4, the five-level delivery-mode-by-birth-order model |
| `R/05_sensitivity_table2.r` | Table 2, the sensitivity summary table |

`R/stepAIC_results.rds` holds the variable names selected in the first AIC round (two character
vectors only, no participant data), so that individual sections can be run without repeating the
full selection.

### Second revision analyses

`R/revision2/` holds the analyses added in the second revision. Every script sources
`R/revision2/_common/setup.R`, which runs the setup and helper blocks of
`R/02_supplementary_forests.r`, so data preparation, variable lists and standardization are those
of the main analysis. Run each script from its own folder; it writes its tables next to itself.
The result tables are included, so the numbers can be checked without the data.

| Folder | Manuscript item |
|---|---|
| `point1_model_equivalence/` | Fig. 3, the final model re-parameterized with one multilevel variable, complete-case and Bayesian |
| `point2_within_stratum_contrasts/` | Supplementary Table 4, effects within birth order strata; the model without the cesarean delivery by first delivery interaction (global effects) |
| `point3_missing_birth_order/` | Supplementary Tables 8 and 9, birth-order missingness model and delta-adjusted tipping-point analysis |
| `point4_delivery_type/` | counts of the four recorded delivery types |
| `point6_cord_glucose_conditioning/` | Supplementary Table 5 and Supplementary Fig. 5, omission models and path analysis |
| `point7_ancestry_heterogeneity/` | Supplementary Fig. 10, ancestry interaction models, likelihood-ratio tests, Cochran's Q |
| `posthoc_power_recalc/`, `power_global_1005/` | Supplementary Table 6, observed power |
| `table2_global_cd_0924/`, `table2_ga_bw_global_1006/` | Table 2, global effects in every sensitivity scenario |
| `supp_v9_1006/`, `supp_v11_table8_1006/`, `supp_v13_table7_1006/` | checks of Supplementary Tables 1, 2, 7 and 8 |
| `fig2a_stars_1006/`, `suppfig_stars_1006/` | Fig. 2a and the supplementary forest plots with the stratum-specific contrasts marked |

Scripts that read a saved JointAI fit need `HAPO_BAYES`, the directory holding the fits
(`Sys.setenv(HAPO_BAYES = "/path/to/fits")`). The fits are not included, for the reason given
below; the `*_fit_jointai_*.R` scripts refit the birth-type and no-interaction joint models
(JointAI 1.1.0, 3 chains, 5,300 iterations, thinning 5, seed 2020).

### Bayesian models

`01_main_pipeline.r` loads two saved JointAI fits (`imputalt.rData`, `imputalt_osszevont.rData`).
These are **not included**: a JointAI object stores the imputed participant-level data, which is
controlled-access. Set `REGENERATE_BAYES <- TRUE` in the script to refit them from the source
data. Refitting takes roughly one hour per model (3 chains, 5,300 iterations, thinning 5,
seed 2020) and requires JAGS.

## Expected output

`outputs/` contains the coefficient tables the scripts produce, so that the analysis can be
checked without access to the data:

- `all_completecase_coefs.csv` — every complete-case model in the paper, one row per coefficient
- `sensitivity_summary_table.csv` — the Table 2 source values
- `supp_fig4_5level_estimates.csv` — the five-level model of Supplementary Fig. 4

## Environment

R 4.2.3 (second revision JointAI fits: R 4.5.2). Packages: `MASS`, `broom`, `dplyr`, `grid`, `pROC`, `ROCR`, `forestploter`, `tableone`,
`flextable`, `officer`, `performance`, `DHARMa`, `pwr`, `sensitivity`, and `JointAI` (with JAGS)
for the Bayesian sections only.

## Citation

See `CITATION.cff`.

## License

MIT, see `LICENSE`.
