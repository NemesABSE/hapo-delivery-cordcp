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

R 4.2.3. Packages: `MASS`, `broom`, `dplyr`, `grid`, `pROC`, `ROCR`, `forestploter`, `tableone`,
`flextable`, `officer`, `performance`, `DHARMa`, `pwr`, `sensitivity`, and `JointAI` (with JAGS)
for the Bayesian sections only.

## Citation

See `CITATION.cff`.

## License

MIT, see `LICENSE`.
