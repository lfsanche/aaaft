# aaaft: rank-based IV inference for accelerated failure time models in Stata

[![Stata 14+](https://img.shields.io/badge/Stata-14%2B-1a5f91)](https://www.stata.com)
[![Version](https://img.shields.io/badge/version-1.1.0-2e7d32)](#version-history)
[![License: MIT](https://img.shields.io/badge/license-MIT-lightgrey)](LICENSE)

`aaaft` estimates the effect of an endogenous treatment on a duration. It works with accelerated failure time (AFT) models and a binary instrument. Its confidence sets are exact and remain valid when the instrument is weak. It does not assume a distribution for the baseline duration or for the frailty.

The command implements the rank inference of Acharya, Khalaf, Menon, Voia, Yazbeck, and Wensley (2026), Canadian Journal of Economics. That method extends the distribution free IV approach of Andrews and Marmer (2008) to durations with multi level unobserved heterogeneity.

Keywords: Stata, econometrics, instrumental variables, weak instruments, duration models, survival analysis, accelerated failure time, rank tests, Hodges Lehmann Sen estimator, health economics.

## Authors

| Author | Affiliation | Contact |
|---|---|---|
| Leonardo Sánchez-Aragón | Facultad de Ciencias Sociales y Humanísticas, Escuela Superior Politécnica del Litoral (ESPOL), Guayaquil, Ecuador | [lfsanche@espol.edu.ec](mailto:lfsanche@espol.edu.ec) |
| Anand Acharya (deceased) | Carleton University, Ottawa, Canada | |
| Kusum Menon | University of Ottawa, Ottawa, Canada | [KMenon@cheo.on.ca](mailto:KMenon@cheo.on.ca) |
| Marcel Voia | Université d'Orléans, France, and University of Bucharest, Romania | [marcel.voia@univ-orleans.fr](mailto:marcel.voia@univ-orleans.fr) |

Questions about the code and bug reports are welcome as [GitHub issues](https://github.com/lfsanche/aaaft/issues) or by email to Leonardo Sánchez-Aragón.

## What the command reports

For each score, `aaaft` reports two rank results and, as a comparator, two least squares results.

1. The Hodges Lehmann Sen (HLS) point estimate, the grid value that minimizes the rank statistic.
2. The rank confidence set. Its critical value is simulated from random permutations, so the set is exact and robust to weak instruments.
3. The least squares IV point estimate.
4. The least squares confidence set from the generalized Anderson Rubin statistic. It is exact only if the errors follow the distribution of the score.

The scores are the logistic, the Gumbel, and the exponential. With the option `fo()`, the command adds an empirical quantile score built from draws of the error that the user supplies. Without an instrument, the command reports a rank only benchmark that ignores endogeneity, for the same scores, including the empirical quantile score with `fo()`. The gap between the two estimates shows the direction and size of the endogeneity bias.

The command treats every duration as observed and does not adjust for censoring. It suits data where all durations are complete, as in the ICU application of the paper. With censored durations the rank test keeps its size, but the point estimate is biased.

## Installation

The command needs Stata 14 or later and the `moremata` package.

```stata
ssc install moremata
net install aaaft, from("https://raw.githubusercontent.com/lfsanche/aaaft/main/") replace
```

To also download the example data and do file to the current folder:

```stata
net get aaaft, from("https://raw.githubusercontent.com/lfsanche/aaaft/main/")
```

To install by hand, copy `aaaft.ado`, `aaaft_export.ado`, and `aaaft.sthlp` to your PERSONAL ado folder. Type `sysdir` in Stata to see where it is. The two ado files must sit in the same folder.

Check the installation with

```stata
which aaaft
help aaaft
```

## Quick start

The file `data/aaaft_example.dta` holds simulated data with the same variable names as the ICU application. The true effect of the severity score on the log length of stay is 0.20. The severity score is endogenous, because frail patients receive lower scores and stay longer.

```stata
use "data/aaaft_example.dta", clear

* Rank IV estimates, with reproducible confidence sets
aaaft lnlosh PRISM hsp1-hsp6 Agecat Cardiac Nchrnx PrevICU, ///
    instrument(Trauma) gridmin(0.01) gridmax(0.40) gridstep(0.005) ///
    reps(1000) seed(12345)

* Rank only benchmark, which ignores endogeneity
aaaft lnlosh PRISM hsp1-hsp6 Agecat Cardiac Nchrnx PrevICU, ///
    gridmin(0.01) gridmax(0.40) gridstep(0.005)
```

The IV estimate should be close to the true value of 0.20, and the rank only benchmark close to 0.08. The do file `examples/aaaft_example.do` runs these steps and exports the results to a CSV file.

## Syntax

```stata
aaaft depvar treatvar covariates [if] [in] [, options]
```

| Option | Description |
|---|---|
| `instrument(varname)` | Binary instrument. Without it, the command reports the rank only benchmark. |
| `gridmin(#)`, `gridmax(#)`, `gridstep(#)` | Grid of candidate values for the coefficient. Defaults 0.01, 0.99, and 0.01. |
| `reps(#)` | Simulation draws for the critical values. Default 1000; 5000 or more for published results. |
| `level(#)` | Confidence level. Default 95. |
| `seed(#)` | Seed for the simulations, so the confidence sets can be reproduced. |
| `fo(filename)` | Stata dataset with draws of the error for the empirical quantile score. |
| `fovar(name)` | Variable in the `fo()` file that holds the draws. Default `eps`. |
| `zcenter` | Use the instrument centered at its mean instead of residualized on the covariates. |
| `norank`, `nols` | Skip the rank results or the least squares results. |
| `export(filename)` | Write the results to a CSV file, one row per score. |

The command stores its results in `r()` and `e()`, so they work with `estimates store`, `estimates table`, `esttab`, and `outreg2`. The vector `e(b)` holds the point estimates only. The bounds are stored as scalars, such as `e(CIR_lo_logistic)` and `e(CIR_hi_logistic)`. Type `help aaaft` for the full list.

## Replication of the paper

The folder `replication` reproduces the empirical results of the Stata Journal article for the pediatric ICU data.

1. `calibrate_fo.do` fits a Weibull model with Gamma frailty in each hospital and saves 100,000 draws of the error, `fo_PRISM.dta` and `fo_PIM.dta`.
2. `running_example.do` runs the IV and rank only estimates for PRISM and PIM and exports them.
3. `robust_fo.do` repeats the empirical quantile estimates under other calibrations of the error.
4. `aaaft_replication.ipynb` runs all the steps in one place, compares the results with Table 2 of Acharya et al. (2026), and runs the validation checks of the command. It uses Stata through `pystata`.

Set the paths at the top of each file before running them.

Data availability. The ICU data cover 10,044 admissions to six Canadian pediatric referral hospitals. They were collected with approval from the institutional review board of each participating center, and they cannot be shared in this repository. Researchers who want to replicate the results can contact the authors. The simulated data in `data` let anyone run every command in the repository.

## Monte Carlo simulations

The folder `simulations` holds the notebook behind Section 5 of the article, `aaaft_montecarlo_v2.ipynb`, with its stored results and figures. It uses a Python implementation of the same algorithm, so it runs without Stata. It checks size and bias across instrument strength, endogeneity, and sample size.

## Repository structure

```
aaaft/
  aaaft.ado              the command
  aaaft_export.ado       CSV export, called by export() or run by itself
  aaaft.sthlp            help file
  aaaft.pkg, stata.toc   files for net install
  data/
    aaaft_example.dta    simulated example data
  examples/
    aaaft_example.do     quick start on the simulated data
  replication/           do files and notebook for the ICU application
  simulations/           Monte Carlo notebook, results, and figures
  CITATION.cff
  LICENSE
```

## Variables in the example data

| Variable | Description |
|---|---|
| `lnlosh` | Log length of stay |
| `PRISM` | Severity score, the endogenous treatment |
| `Trauma` | Trauma admission, the instrument |
| `hsp1` to `hsp6` | Hospital indicators |
| `Agecat` | Age category, 1 to 4 |
| `Cardiac` | Cardiac diagnosis |
| `Nchrnx` | Number of chronic conditions |
| `PrevICU` | Previous ICU admission |

## Citation

If you use `aaaft`, please cite the article and the method.

```bibtex
@article{sanchezaragon2026aaaft,
  author  = {S{\'a}nchez-Arag{\'o}n, Leonardo and Acharya, Anand and Menon, Kusum and Voia, Marcel},
  title   = {aaaft: Rank-based instrumental variables inference for accelerated failure time models in {Stata}},
  year    = {2026},
  note    = {Working paper}
}

@article{acharya2026,
  author  = {Acharya, Anand and Khalaf, Lynda and Menon, Kusum and Voia, Marcel and Yazbeck, Myra and Wensley, David},
  title   = {Reducing Bias in Pediatric Critical Care Allocation: Evidence from Canadian Hospitals using Robust Instrumental Variable Duration Models},
  journal = {Canadian Journal of Economics},
  year    = {2026},
  note    = {Forthcoming}
}
```

## Version history

1.1.0 (October 2026). The empirical quantile score from `fo()` is also available in the rank only benchmark, without an instrument.

1.0.0 (October 2026). First public release. Rank IV estimates and simulated rank confidence sets under the logistic, Gumbel, exponential, and empirical quantile scores; least squares IV comparator; rank only benchmark; CSV export.

## License

The code is released under the [MIT License](LICENSE).

## Acknowledgments

We thank Kim Huynh, Steven Lehrer, Pravin Trivedi, and Tiemen Woutersen for their contributions to the development of this project, and the Canadian Critical Care Trials Group for data collection support.
