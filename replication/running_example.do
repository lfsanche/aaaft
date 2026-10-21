* running_example.do
* Replication example for the aaaft command (Stata Journal).
* Estimates the effect of illness severity on ICU length of stay for PRISM
* and PIM, with and without the trauma instrument, and exports the IV results.
* Requires aaaft.ado and aaaft_export.ado (version 1.0.0 or later).
* Run calibrate_fo.do first: it creates fo_PRISM.dta and fo_PIM.dta,
* the error draws for the empirical quantile score.

clear all
version 14

* Dependency, install once:
* ssc install moremata

* Output folder and data. Change these two paths to your own.
cd "path/to/output"
use "path/to/HE_May2023_V13", clear   // ICU data, available from the authors

* Common settings. With reps(5000) most confidence bounds move at most one
* grid step across seeds; the seed makes them reproducible.
local covs  hsp1-hsp6 Agecat Cardiac Nchrnx PrevICU
local grid  gridmin(0.01) gridmax(0.60) gridstep(0.005)
local sim   reps(5000) seed(12345)

* -------------------------------------------------------------------
* PRISM
* -------------------------------------------------------------------
* IV estimates with the empirical quantile score, exported to PRISM_iv.csv
aaaft lnlosh PRISM `covs', instrument(Trauma) `grid' `sim' ///
    fo("fo_PRISM.dta") export("PRISM_iv.csv")
estimates store prism_iv
scalar prism_iv = e(HLS_logistic)

* Rank-only benchmark, no instrument
aaaft lnlosh PRISM `covs', `grid'
estimates store prism_noiv
scalar prism_noiv = e(HLS_logistic)

* -------------------------------------------------------------------
* PIM
* -------------------------------------------------------------------
* IV estimates with the empirical quantile score, exported to PIM_iv.csv
aaaft lnlosh PIM `covs', instrument(Trauma) `grid' `sim' ///
    fo("fo_PIM.dta") export("PIM_iv.csv")
estimates store pim_iv
scalar pim_iv = e(HLS_logistic)

* Rank-only benchmark, no instrument
aaaft lnlosh PIM `covs', `grid'
estimates store pim_noiv
scalar pim_noiv = e(HLS_logistic)

* -------------------------------------------------------------------
* Point estimates side by side
* -------------------------------------------------------------------
estimates table prism_iv prism_noiv pim_iv pim_noiv, b(%6.3f)

* -------------------------------------------------------------------
* Endogeneity gap, logistic score
* -------------------------------------------------------------------
scalar prism_gap = prism_iv - prism_noiv
scalar pim_gap   = pim_iv   - pim_noiv
di as txt "PRISM  IV = " as res %6.3f prism_iv ///
    as txt "   no-IV = " as res %6.3f prism_noiv ///
    as txt "   gap = " as res %6.3f prism_gap
di as txt "PIM    IV = " as res %6.3f pim_iv ///
    as txt "   no-IV = " as res %6.3f pim_noiv ///
    as txt "   gap = " as res %6.3f pim_gap

* -------------------------------------------------------------------
* Read back the exported results
* -------------------------------------------------------------------
* export() calls aaaft_export. It can also be run by itself after an
* estimation, for example: aaaft_export using PRISM_iv.csv, replace
preserve
import delimited using "PRISM_iv.csv", clear
list score hls rank_lo rank_hi ls_point ls_lo ls_hi, noobs
restore
