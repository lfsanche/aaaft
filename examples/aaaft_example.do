* aaaft_example.do
* Quick start for aaaft with simulated data.
* The true effect of PRISM on lnlosh is 0.20. PRISM is endogenous: frail
* patients get lower scores and stay longer, so ignoring the endogeneity
* biases the estimate down. Trauma is the instrument.
* Run from the root of the repository, or after net get from the same folder.

clear all
version 14

capture use "data/aaaft_example.dta", clear
if _rc use "aaaft_example.dta", clear

local covs hsp1-hsp6 Agecat Cardiac Nchrnx PrevICU
local grid gridmin(0.01) gridmax(0.40) gridstep(0.005)

* Rank IV estimates, exported to a CSV file
aaaft lnlosh PRISM `covs', instrument(Trauma) `grid' reps(1000) seed(12345) ///
    export("aaaft_example_results.csv")
estimates store iv

* Rank only benchmark, ignoring the endogeneity of PRISM
aaaft lnlosh PRISM `covs', `grid'
estimates store noiv

estimates table iv noiv, b(%6.3f)
