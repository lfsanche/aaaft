* calibrate_fo.do
* Construye las realizaciones del error para el score de cuantil empírico de aaivaft.
* Modelo: AFT Weibull con fragilidad Gamma, estimado por separado en cada
* hospital (streg no admite ancillary() junto con frailty()).
* Error de ln(T): eps = sigma_g * (W - ln(alpha))
*   W     : valor extremo de mínimo estándar (error de ln T en un Weibull)
*   alpha : fragilidad Gamma con media 1 y varianza theta_g
*   sigma_g = 1/p_g, escala del hospital g
* Con theta = 1 y una sola escala, eps es exactamente logístico.
* Se calibra por separado para PRISM y PIM, ignorando la endogeneidad:
* solo se busca la forma del error, no el coeficiente.

clear all
version 14

local data "path/to/HE_May2023_V13"   // ICU data, available from the authors
local out  "path/to/output"
local B    100000

foreach sev in PRISM PIM {

    * -------------------------------------------------------------
    * 1. Estimación del modelo Weibull con fragilidad Gamma
    * -------------------------------------------------------------
    use "`data'", clear

    * Las seis dummies de hospital deben sumar uno
    egen byte hs_total = rowtotal(hsp1-hsp6)
    assert hs_total == 1

    * Todas las estancias se observan completas (sin censura)
    gen double t = exp(lnlosh)
    stset t

    * Un modelo por hospital: escala sigma_g y varianza de fragilidad theta_g
    quietly count
    local ntot = r(N)
    forvalues g = 1/6 {
        quietly streg `sev' Agecat Cardiac Nchrnx PrevICU if hsp`g' == 1, ///
            distribution(weibull) frailty(gamma) time
        scalar sig`g'   = exp(-_b[/ln_p])
        scalar theta`g' = e(theta)
        scalar sh`g'    = e(N) / `ntot'
    }

    di as txt _n "Calibration for `sev'"
    forvalues g = 1/6 {
        di as txt "  hospital `g': sigma = " as res %7.4f sig`g' ///
            as txt "   theta = " as res %7.4f theta`g' ///
            as txt "   share = " as res %6.3f sh`g'
    }

    * -------------------------------------------------------------
    * 2. Simulación de B realizaciones del error
    * -------------------------------------------------------------
    clear
    set obs `B'
    set seed 2026

    * Hospital de cada realización, según la participación en la muestra
    gen double u = runiform()
    gen byte g = 1
    local cum = 0
    forvalues k = 1/5 {
        local cum = `cum' + sh`k'
        quietly replace g = `k' + 1 if u > `cum'
    }

    gen double sig = .
    gen double th  = .
    forvalues k = 1/6 {
        quietly replace sig = sig`k'   if g == `k'
        quietly replace th  = theta`k' if g == `k'
    }

    gen double W     = ln(-ln(runiform()))
    gen double alpha = cond(th < 1e-6, 1, rgamma(1/th, th))
    gen double eps   = sig * (W - ln(alpha))

    label variable eps "Simulated AFT error, Weibull-gamma frailty, `sev'"
    keep eps
    summarize eps, detail
    save "`out'\fo_`sev'.dta", replace
}

* -------------------------------------------------------------
* 3. aaivaft con el score de cuantil empírico
* -------------------------------------------------------------
use "`data'", clear

aaivaft lnlosh PRISM hsp1-hsp6 Agecat Cardiac Nchrnx PrevICU, ///
    instrument(Trauma) gridmin(0.01) gridmax(0.60) gridstep(0.005) ///
    reps(5000) seed(12345) nols fo("`out'\fo_PRISM.dta")

aaivaft lnlosh PIM hsp1-hsp6 Agecat Cardiac Nchrnx PrevICU, ///
    instrument(Trauma) gridmin(0.01) gridmax(0.60) gridstep(0.005) ///
    reps(5000) seed(12345) nols fo("`out'\fo_PIM.dta")
