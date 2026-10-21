* robust_fo.do
* Robustez del score de cuantil empírico de aaaft a la calibración de Fo.
* Requiere aaaft versión 1.0.0 o posterior.
* Calibraciones (error de ln T: eps = sigma_g * (W - ln(alpha))):
*   base    : sigma_g y theta_g estimados por hospital (como calibrate_fo.do)
*   pooled  : un solo modelo, sigma y theta comunes a todos los hospitales
*   theta1  : sigma_g por hospital, theta = 1 (mezcla de logísticas)
*   theta_x2: sigma_g por hospital, theta_g al doble
*   theta_x05: sigma_g por hospital, theta_g a la mitad
*   split   : base calibrado en una mitad de la muestra, aaaft en la otra
*             (conserva la exactitud; se compara con el logístico en esa mitad)

clear all
version 14

local data "path/to/HE_May2023_V13"   // ICU data, available from the authors
local out  "path/to/output"
global B   100000
local covs Agecat Cardiac Nchrnx PrevICU
local opts instrument(Trauma) gridmin(0.01) gridmax(0.60) gridstep(0.005) reps(5000) seed(12345) nols

* -------------------------------------------------------------
* Programa: calibra sigma_g, theta_g y pesos por hospital
* Deja escalares sig1-sig6, th1-th6, sh1-sh6
* -------------------------------------------------------------
capture program drop calib_hosp
program define calib_hosp
    args sev covs cond
    quietly count `cond'
    local ntot = r(N)
    forvalues g = 1/6 {
        quietly streg `sev' `covs' `cond' & hsp`g' == 1, ///
            distribution(weibull) frailty(gamma) time
        scalar sig`g' = exp(-_b[/ln_p])
        scalar th`g'  = e(theta)
        scalar sh`g'  = e(N) / `ntot'
    }
end

* -------------------------------------------------------------
* Programa: simula $B realizaciones de eps y las guarda
* Usa los escalares sig1-sig6, th1-th6, sh1-sh6
* -------------------------------------------------------------
capture program drop sim_fo
program define sim_fo
    args file
    preserve
    clear
    quietly set obs $B
    set seed 2026
    quietly gen double u = runiform()
    quietly gen byte g = 1
    local cum = 0
    forvalues k = 1/5 {
        local cum = `cum' + sh`k'
        quietly replace g = `k' + 1 if u > `cum'
    }
    quietly gen double sig = .
    quietly gen double th  = .
    forvalues k = 1/6 {
        quietly replace sig = sig`k' if g == `k'
        quietly replace th  = th`k'  if g == `k'
    }
    quietly gen double W     = ln(-ln(runiform()))
    quietly gen double alpha = cond(th < 1e-6, 1, rgamma(1/th, th))
    quietly gen double eps   = sig * (W - ln(alpha))
    keep eps
    quietly save `"`file'"', replace
    restore
end

* -------------------------------------------------------------
* Resultados
* -------------------------------------------------------------
tempfile res
postfile h str6 sev str10 calib hls_q lo_q hi_q hls_l lo_l hi_l using `res'

foreach sev in PRISM PIM {

    use "`data'", clear
    gen double t = exp(lnlosh)
    quietly stset t

    * --- base: por hospital
    calib_hosp `sev' "`covs'" "if 1"
    forvalues g = 1/6 {
        scalar bsig`g' = sig`g'
        scalar bth`g'  = th`g'
    }
    sim_fo "`out'\rob_`sev'_base.dta"

    * --- theta = 1
    forvalues g = 1/6 {
        scalar th`g' = 1
    }
    sim_fo "`out'\rob_`sev'_theta1.dta"

    * --- theta al doble y a la mitad
    forvalues g = 1/6 {
        scalar th`g' = 2 * bth`g'
    }
    sim_fo "`out'\rob_`sev'_thx2.dta"
    forvalues g = 1/6 {
        scalar th`g' = 0.5 * bth`g'
    }
    sim_fo "`out'\rob_`sev'_thx05.dta"

    * --- pooled: un solo modelo con dummies de hospital
    quietly streg `sev' `covs' hsp2-hsp6, ///
        distribution(weibull) frailty(gamma) time
    forvalues g = 1/6 {
        scalar sig`g' = exp(-_b[/ln_p])
        scalar th`g'  = e(theta)
    }
    di as txt "`sev' pooled: sigma = " as res %6.4f sig1 as txt "  theta = " as res %6.4f th1
    sim_fo "`out'\rob_`sev'_pooled.dta"

    * --- aaaft con cada calibración en la muestra completa
    foreach c in base pooled theta1 thx2 thx05 {
        quietly aaaft lnlosh `sev' hsp1-hsp6 `covs', `opts' ///
            fo("`out'\rob_`sev'_`c'.dta")
        post h ("`sev'") ("`c'") (e(HLS_empq)) (e(CIR_lo_empq)) (e(CIR_hi_empq)) ///
            (e(HLS_logistic)) (e(CIR_lo_logistic)) (e(CIR_hi_logistic))
    }

    * --- split: calibrar en una mitad, estimar en la otra
    set seed 777
    gen byte half = runiform() < 0.5
    calib_hosp `sev' "`covs'" "if half == 1"
    sim_fo "`out'\rob_`sev'_split.dta"
    quietly aaaft lnlosh `sev' hsp1-hsp6 `covs' if half == 0, `opts' ///
        fo("`out'\rob_`sev'_split.dta")
    post h ("`sev'") ("split") (e(HLS_empq)) (e(CIR_lo_empq)) (e(CIR_hi_empq)) ///
        (e(HLS_logistic)) (e(CIR_lo_logistic)) (e(CIR_hi_logistic))
}
postclose h

use `res', clear
format hls_q lo_q hi_q hls_l lo_l hi_l %6.3f
list, noobs sepby(sev)
