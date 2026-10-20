*! version 1.0.0  19oct2026
*! aaivaft: rank-based IV inference for accelerated failure time models
*! Authors: Sanchez-Aragon, Acharya, Menon, and Voia
*! Based on: Acharya, Khalaf, Menon, Voia, Yazbeck, and Wensley (2026), Canadian Journal of Economics, forthcoming
*! Requires: moremata (ssc install moremata) and aaivaft_export.ado

program define aaivaft, rclass
	version 14.0

	syntax varlist(min=2 numeric) 		///
		[if] [in]						///
		[, 								///
		INSTRument(varname numeric)		///
		GRIDmin(real 0.01)				///
		GRIDmax(real 0.99)				///
		GRIDstep(real 0.01)				///
		REPS(integer 1000)				///
		Level(cilevel)					///
		noRANK							///
		noLS							///
        EXPORT(string)                  ///
		SEED(string)					///
		FO(string)						///
		FOVar(name)						///
		ZCENTer							///
		]

	marksample touse

	// Seed for the rank and LS simulations
	if `"`seed'"' != "" {
		set seed `seed'
	}

	// Empirical quantile score from user-supplied error draws
	local has_fo 0
	if `"`fo'"' != "" {
		if "`fovar'" == "" local fovar eps
		capture confirm file `"`fo'"'
		if _rc {
			capture confirm file `"`fo'.dta"'
			if _rc {
				di as error `"fo(): file `fo' not found."'
				exit 601
			}
			local fo `"`fo'.dta"'
		}
		preserve
		quietly use `"`fo'"', clear
		capture confirm numeric variable `fovar'
		if _rc {
			di as error "fo(): numeric variable `fovar' not found in `fo'."
			di as error "Use fovar() to name the variable that holds the draws."
			exit 111
		}
		quietly count if !missing(`fovar')
		if r(N) < 1000 {
			di as error "fo(): at least 1,000 nonmissing draws are needed."
			exit 2001
		}
		local fo_n = r(N)
		mata: _aaivaft_fodraws = st_data(., "`fovar'")
		mata: _aaivaft_fodraws = select(_aaivaft_fodraws, _aaivaft_fodraws :< .)
		restore
		local has_fo 1
	}
	else if "`fovar'" != "" {
		di as error "fovar() requires fo()."
		exit 198
	}

	// -------------------------------------------------------------------------
	// Parse varlist: depvar treatvar [covariates]
	// -------------------------------------------------------------------------
	local depvar   : word 1 of `varlist'
	local treatvar : word 2 of `varlist'
	local total_vars : word count `varlist'

	if `total_vars' < 3 {
		di as error "You must supply at least one covariate after the treatment variable."
		exit 102
	}

	local covvars ""
	forvalues i = 3/`total_vars' {
		local v : word `i' of `varlist'
		local covvars `covvars' `v'
	}

	// If no instrument: force noLS (LS CI requires IV)
	local has_inst = ("`instrument'" != "")
	if !`has_inst' {
		local ls "nols"
	}

	// -------------------------------------------------------------------------
	// Dependency check
	// -------------------------------------------------------------------------
	capture mata: mm_ranks((1\2\3), 1, 0)
	if _rc {
		di as error "This command requires the -moremata- package."
		di as error "Please install it: {stata ssc install moremata}"
		exit 499
	}

	// -------------------------------------------------------------------------
	// Input validation
	// -------------------------------------------------------------------------
	if `gridmin' >= `gridmax' {
		di as error "gridmin() must be strictly less than gridmax()."
		exit 198
	}
	if `gridstep' <= 0 {
		di as error "gridstep() must be positive."
		exit 198
	}
	if `reps' < 100 {
		di as error "reps() must be at least 100."
		exit 198
	}

	// -------------------------------------------------------------------------
	// Count observations in sample
	// -------------------------------------------------------------------------
	quietly count if `touse'
	local nobs = r(N)
	if `nobs' < 10 {
		di as error "Fewer than 10 observations in sample. Cannot proceed."
		exit 2001
	}

	// -------------------------------------------------------------------------
	// Drop covariates that are constant within the sample (fix for if/in)
	// -------------------------------------------------------------------------
	local covvars_clean ""
	local dropped ""
	foreach v of local covvars {
		quietly summarize `v' if `touse'
		if r(sd) > 0 {
			local covvars_clean `covvars_clean' `v'
		}
		else {
			local dropped `dropped' `v'
		}
	}
	if "`dropped'" != "" {
		di as txt " Note: dropped constant covariates in this sample: " ///
			as res "`dropped'"
	}
	if "`covvars_clean'" == "" {
		di as error "All covariates are constant in the selected sample."
		exit 198
	}
	local covvars `covvars_clean'

	// Treatment and instrument must vary in the estimation sample
	quietly summarize `treatvar' if `touse'
	if r(sd) == 0 | r(sd) == . {
		di as error "Treatment variable `treatvar' is constant in the selected sample."
		exit 459
	}
	if `has_inst' {
		quietly summarize `instrument' if `touse'
		if r(sd) == 0 | r(sd) == . {
			di as error "Instrument `instrument' is constant in the selected sample."
			exit 459
		}
	}

	// -------------------------------------------------------------------------
	// Header
	// -------------------------------------------------------------------------
	di _n as txt "{hline 60}"
	di as txt " Rank IV inference for AFT duration models (aaivaft)"
	di as txt "{hline 60}"
	di as txt " Dependent variable  : " as res "`depvar'"
	di as txt " Treatment variable  : " as res "`treatvar'"
	if `has_inst' {
		di as txt " Instrument          : " as res "`instrument'"
	}
	else {
		di as txt " Instrument          : " as res "(none — rank estimates only)"
	}
	di as txt " Covariates          : " as res "`covvars'"
	di as txt " Observations        : " as res "`nobs'"
	di as txt " Beta grid           : " as res "`gridmin'" as txt " to " ///
		as res "`gridmax'" as txt " (step " as res "`gridstep'" as txt ")"
	if `has_inst' {
		di as txt " Simulations (reps)  : " as res "`reps'"
		di as txt " Confidence level    : " as res "`level'" as txt "%"
	}
	di as txt "{hline 60}"

	// -------------------------------------------------------------------------
	// Build matrices for Mata
	// -------------------------------------------------------------------------
	tempname _aaivaft_y _aaivaft_Y _aaivaft_X

	quietly {
		mkmat `depvar'   if `touse', matrix(`_aaivaft_y')
		mkmat `treatvar' if `touse', matrix(`_aaivaft_Y')

		local first_cov : word 1 of `covvars'
		mkmat `first_cov' if `touse', matrix(`_aaivaft_X')
		local ncov : word count `covvars'
		if `ncov' > 1 {
			forvalues i = 2/`ncov' {
				local cv : word `i' of `covvars'
				tempname _tmp
				mkmat `cv' if `touse', matrix(`_tmp')
				matrix `_aaivaft_X' = `_aaivaft_X', `_tmp'
			}
		}
	}

	// Initialize HLS scalars (set missing in case norank is specified)
	tempname _hls_bl _hls_bg _hls_be _hls_bq
	scalar `_hls_bl' = .
	scalar `_hls_bg' = .
	scalar `_hls_be' = .
	scalar `_hls_bq' = .

	// Rank confidence sets (lo, hi, edge flag, gap flag, critical value)
	tempname _cr_lo_l _cr_hi_l _cr_ed_l _cr_gp_l _cr_gc_l
	tempname _cr_lo_g _cr_hi_g _cr_ed_g _cr_gp_g _cr_gc_g
	tempname _cr_lo_e _cr_hi_e _cr_ed_e _cr_gp_e _cr_gc_e
	tempname _cr_lo_q _cr_hi_q _cr_ed_q _cr_gp_q _cr_gc_q
	foreach s in l g e q {
		foreach q in lo hi ed gp gc {
			scalar `_cr_`q'_`s'' = .
		}
	}

	// Build instrument matrix once — used by both rank and LS blocks
	if `has_inst' {
		tempname _aaivaft_z
		quietly mkmat `instrument' if `touse', matrix(`_aaivaft_z')
	}

	scalar _aaivaft_gmin  = `gridmin'
	scalar _aaivaft_gmax  = `gridmax'
	scalar _aaivaft_gstep = `gridstep'
	scalar _aaivaft_reps  = `reps'
	scalar _aaivaft_level = `level'

	// -------------------------------------------------------------------------
	// Rank estimator
	// -------------------------------------------------------------------------
	// Fo() only applies to the rank IV estimates
	if `has_fo' & (!`has_inst' | "`rank'" != "") {
		di as txt " Note: fo() is used only with instrument() and without norank; ignored."
		local has_fo 0
	}
	local sq ""
	if `has_fo' local sq "q"

	if "`rank'" == "" {
		di as txt _n " Computing rank estimates..."

		if `has_inst' {
			mata: _aaivaft_rankinf(		///
				"`_aaivaft_y'",			///
				"`_aaivaft_Y'",			///
				"`_aaivaft_z'",			///
				"`_aaivaft_X'",			///
				"_aaivaft_gmin",			///
				"_aaivaft_gmax",			///
				"_aaivaft_gstep",			///
				"_aaivaft_reps",			///
				"_aaivaft_level",			///
				`has_fo',				///
				`=("`zcenter'" == "")'	///
			)
		}
		else {
			// No-IV version: uses treatment directly (OLS-rank)
			mata: _aaivaft_rankinf_noinst(	///
				"`_aaivaft_y'",				///
				"`_aaivaft_Y'",				///
				"`_aaivaft_X'",				///
				"_aaivaft_gmin",				///
				"_aaivaft_gmax",				///
				"_aaivaft_gstep"				///
			)
		}
	}

	// Save HLS results before LS overwrites r()
	if "`rank'" == "" {
		scalar `_hls_bl' = r(beta_logistic)
		scalar `_hls_bg' = r(beta_gumbel)
		scalar `_hls_be' = r(beta_exponential)
		if `has_fo' scalar `_hls_bq' = r(beta_empq)
		if `has_inst' {
			foreach s in l g e `sq' {
				if "`s'" == "l" local sn logistic
				if "`s'" == "g" local sn gumbel
				if "`s'" == "e" local sn exponential
				if "`s'" == "q" local sn empq
				scalar `_cr_lo_`s'' = r(cir_l_`sn')
				scalar `_cr_hi_`s'' = r(cir_u_`sn')
				scalar `_cr_ed_`s'' = r(cir_edge_`sn')
				scalar `_cr_gp_`s'' = r(cir_gap_`sn')
				scalar `_cr_gc_`s'' = r(gcrit_`sn')
			}
		}
	}

	// -------------------------------------------------------------------------
	// LS confidence sets (only available with instrument)
	// -------------------------------------------------------------------------
	if `has_inst' & "`ls'" == "" {
		di as txt _n " Computing least-squares IV confidence sets..."
		mata: _aaivaft_leastsquares(		///
			"`_aaivaft_y'",				///
			"`_aaivaft_Y'",				///
			"`_aaivaft_z'",				///
			"`_aaivaft_X'",				///
			"_aaivaft_reps",				///
			"_aaivaft_level"				///
		)
	}

	// Save LS CI results (may be missing if noLS or no instrument)
	tempname _ci_l_l _ci_u_l _ci_l_g _ci_u_g _ci_l_e _ci_u_e _ls_b
	scalar `_ls_b' = .
	if `has_inst' & "`ls'" == "" {
		scalar `_ls_b'   = r(beta_ls)
		scalar `_ci_l_l' = r(ci_l_logistic)
		scalar `_ci_u_l' = r(ci_u_logistic)
		scalar `_ci_l_g' = r(ci_l_gumbel)
		scalar `_ci_u_g' = r(ci_u_gumbel)
		scalar `_ci_l_e' = r(ci_l_exponential)
		scalar `_ci_u_e' = r(ci_u_exponential)
	}
	else {
		scalar `_ci_l_l' = .
		scalar `_ci_u_l' = .
		scalar `_ci_l_g' = .
		scalar `_ci_u_g' = .
		scalar `_ci_l_e' = .
		scalar `_ci_u_e' = .
	}

	// Restore HLS results into r() after LS may have overwritten them
	if "`rank'" == "" {
		return scalar beta_logistic  = `_hls_bl'
		return scalar beta_gumbel    = `_hls_bg'
		return scalar beta_exponential = `_hls_be'
		if `has_fo' {
			return scalar beta_empq = `_hls_bq'
			return scalar fo_draws  = `fo_n'
			return local  fo `"`fo'"'
		}
		if `has_inst' {
			foreach s in l g e `sq' {
				if "`s'" == "l" local sn logistic
				if "`s'" == "g" local sn gumbel
				if "`s'" == "e" local sn exponential
				if "`s'" == "q" local sn empq
				return scalar cir_l_`sn'    = `_cr_lo_`s''
				return scalar cir_u_`sn'    = `_cr_hi_`s''
				return scalar cir_edge_`sn' = `_cr_ed_`s''
				return scalar cir_gap_`sn'  = `_cr_gp_`s''
				return scalar gcrit_`sn'    = `_cr_gc_`s''
			}
		}
	}

	if `has_inst' & "`ls'" == "" {
		return scalar beta_ls = `_ls_b'
		return scalar ci_l_logistic    = `_ci_l_l'
		return scalar ci_u_logistic    = `_ci_u_l'
		return scalar ci_l_gumbel      = `_ci_l_g'
		return scalar ci_u_gumbel      = `_ci_u_g'
		return scalar ci_l_exponential = `_ci_l_e'
		return scalar ci_u_exponential = `_ci_u_e'
	}

	// -------------------------------------------------------------------------
	// Results table
	// -------------------------------------------------------------------------
	di _n as txt "{hline 60}"
	di as txt "  Results"
	di as txt "{hline 60}"

	if "`rank'" == "" {
		if `has_inst' {
			di as txt _n "  Hodges-Lehmann-Sen point estimates (Rank IV):"
		}
		else {
			di as txt _n "  Hodges-Lehmann-Sen point estimates (Rank, no IV):"
		}
		di as txt "  {hline 48}"
		di as txt "  {col 4}Score{col 28}Beta"
		di as txt "  {hline 48}"
		di as txt "  {col 4}Logistic" ///
			"{col 28}" as res %8.4f `_hls_bl'
		di as txt "  {col 4}Gumbel (Extreme Value)" ///
			"{col 28}" as res %8.4f `_hls_bg'
		di as txt "  {col 4}Exponential" ///
			"{col 28}" as res %8.4f `_hls_be'
		if `has_fo' {
			di as txt "  {col 4}Empirical quantile (fo)" ///
				"{col 28}" as res %8.4f `_hls_bq'
		}
		di as txt "  {hline 48}"
	}

	if "`rank'" == "" & `has_inst' {
		di as txt _n "  `level'% Confidence Sets (Rank IV):"
		di as txt "  {hline 48}"
		di as txt "  {col 4}Score{col 20}Lower{col 34}Upper"
		di as txt "  {hline 48}"
		local anyflag 0
		foreach s in l g e `sq' {
			if "`s'" == "l" local lab "Logistic"
			if "`s'" == "g" local lab "Gumbel"
			if "`s'" == "e" local lab "Exponential"
			if "`s'" == "q" local lab "Emp. quantile"
			local fl ""
			if `_cr_ed_`s'' == 1 local fl "`fl'*"
			if `_cr_gp_`s'' == 1 local fl "`fl'+"
			if "`fl'" != "" local anyflag 1
			if missing(`_cr_lo_`s'') {
				di as txt "  {col 4}`lab'{col 20}" as res "(empty on grid)"
			}
			else {
				di as txt "  {col 4}`lab'" ///
					"{col 20}" as res %8.4f `_cr_lo_`s'' ///
					"{col 34}" as res %8.4f `_cr_hi_`s'' ///
					as txt " `fl'"
			}
		}
		di as txt "  {hline 48}"
		if `anyflag' {
			di as txt "  * set reaches the grid edge: widen gridmin()/gridmax()"
			di as txt "  + set is not connected on the grid: bounds are min/max"
		}
	}

	if `has_inst' & "`ls'" == "" {
		di as txt _n "  Least-squares IV point estimate:" ///
			"{col 37}" as res %8.4f `_ls_b'
		di as txt _n "  `level'% Confidence Sets (Least-Squares IV):"
		di as txt "  {hline 48}"
		di as txt "  {col 4}Score{col 20}Lower{col 34}Upper"
		di as txt "  {hline 48}"
		di as txt "  {col 4}Logistic" ///
			"{col 20}" as res %8.4f `_ci_l_l' ///
			"{col 34}" as res %8.4f `_ci_u_l'
		di as txt "  {col 4}Gumbel" ///
			"{col 20}" as res %8.4f `_ci_l_g' ///
			"{col 34}" as res %8.4f `_ci_u_g'
		di as txt "  {col 4}Exponential" ///
			"{col 20}" as res %8.4f `_ci_l_e' ///
			"{col 34}" as res %8.4f `_ci_u_e'
		di as txt "  {hline 48}"
	}

	di as txt "{hline 60}"
	if `has_inst' {
		if "`rank'" == "" {
			di as txt " Note: Rank confidence sets are exact and robust to weak"
			di as txt "       instruments."
		}
		if "`ls'" == "" {
			di as txt " Note: LS confidence sets are exact only if the errors follow"
			di as txt "       the distribution of each score."
		}
		di as txt " Instrument: `instrument'"
		if "`zcenter'" != "" & "`rank'" == "" {
			di as txt " Rank statistic uses the centered instrument (zcenter)."
		}
		if `has_fo' {
			di as txt `" Empirical quantile score: `fo_n' draws of `fovar' from `fo'"'
		}
	}
	else {
		di as txt " Note: No instrument supplied. Only rank point estimates shown."
		di as txt " To obtain confidence sets, supply instrument() option."
	}
	di as txt "{hline 60}"

	scalar drop _aaivaft_gmin _aaivaft_gmax _aaivaft_gstep _aaivaft_reps _aaivaft_level

	// -------------------------------------------------------------------------
	// Delegate ereturn post to eclass helper (avoids r(151) conflict)
	// Values are passed via fixed-name global scalars (prefix _aaivaft_ep_)
	// to avoid the tempname-expansion problem in option parsing.
	// -------------------------------------------------------------------------
	local do_rank = ("`rank'" == "")
	local do_ls   = (`has_inst' & "`ls'" == "")

	if !`do_rank' & !`do_ls' {
		di as txt " Note: no estimates posted to e() (norank and noLS both specified)."
	}
	else {
		// Stage results in fixed-name scalars for the helper to read
		scalar _aaivaft_ep_hls_l  = `_hls_bl'
		scalar _aaivaft_ep_hls_g  = `_hls_bg'
		scalar _aaivaft_ep_hls_e  = `_hls_be'
		scalar _aaivaft_ep_ci_l_l = `_ci_l_l'
		scalar _aaivaft_ep_ci_u_l = `_ci_u_l'
		scalar _aaivaft_ep_ci_l_g = `_ci_l_g'
		scalar _aaivaft_ep_ci_u_g = `_ci_u_g'
		scalar _aaivaft_ep_ci_l_e = `_ci_l_e'
		scalar _aaivaft_ep_ci_u_e = `_ci_u_e'
		scalar _aaivaft_ep_nobs    = `nobs'
		scalar _aaivaft_ep_has_inst = `has_inst'
		scalar _aaivaft_ep_do_ls   = `do_ls'
		scalar _aaivaft_ep_ls_b    = `_ls_b'
		scalar _aaivaft_ep_do_rank = `do_rank'
		scalar _aaivaft_ep_level   = `level'
		foreach s in l g e q {
			scalar _aaivaft_ep_cr_lo_`s' = `_cr_lo_`s''
			scalar _aaivaft_ep_cr_hi_`s' = `_cr_hi_`s''
		}
		scalar _aaivaft_ep_hls_q  = `_hls_bq'
		scalar _aaivaft_ep_has_fo = `has_fo'
		global AAIVAFT_ep_fo `"`fo'"'

		_aaivaft_epost `depvar' `treatvar' `instrument' `touse'
	}

	// Drop the stored draws and the staged file name
	capture mata: mata drop _aaivaft_fodraws
	global AAIVAFT_ep_fo

    // -------------------------------------------------------------------------
    // Optional export to .csv (calls aaivaft_export.ado; reads e())
    // -------------------------------------------------------------------------
    if `"`export'"' != "" {
        if `do_rank' | `do_ls' {
            aaivaft_export using `"`export'"', replace
        }
        else {
            di as txt " Note: nothing posted to e(); export() skipped."
        }
    }

end


// =============================================================================
// eclass helper: posts results to e() — called at end of aaivaft
// Reads values from fixed-name scalars _aaivaft_ep_* set by the caller.
// Positional args: depvar treatvar [instrument] touse
// =============================================================================
program define _aaivaft_epost, eclass
	version 14.0

	args depvar treatvar instrument touse

	// Read staged scalars
	local nobs     = scalar(_aaivaft_ep_nobs)
	local has_inst = scalar(_aaivaft_ep_has_inst)
	local do_ls    = scalar(_aaivaft_ep_do_ls)
	local do_rank  = scalar(_aaivaft_ep_do_rank)
	local lvl      = scalar(_aaivaft_ep_level)

	local hls_l  = scalar(_aaivaft_ep_hls_l)
	local hls_g  = scalar(_aaivaft_ep_hls_g)
	local hls_e  = scalar(_aaivaft_ep_hls_e)
	local ci_l_l = scalar(_aaivaft_ep_ci_l_l)
	local ci_u_l = scalar(_aaivaft_ep_ci_u_l)
	local ci_l_g = scalar(_aaivaft_ep_ci_l_g)
	local ci_u_g = scalar(_aaivaft_ep_ci_u_g)
	local ci_l_e = scalar(_aaivaft_ep_ci_l_e)
	local ci_u_e = scalar(_aaivaft_ep_ci_u_e)
	local ls_b   = scalar(_aaivaft_ep_ls_b)
	scalar drop _aaivaft_ep_ls_b
	local hls_q  = scalar(_aaivaft_ep_hls_q)
	local has_fo = scalar(_aaivaft_ep_has_fo)
	scalar drop _aaivaft_ep_hls_q _aaivaft_ep_has_fo
	foreach s in l g e q {
		local cr_lo_`s' = scalar(_aaivaft_ep_cr_lo_`s')
		local cr_hi_`s' = scalar(_aaivaft_ep_cr_hi_`s')
		scalar drop _aaivaft_ep_cr_lo_`s' _aaivaft_ep_cr_hi_`s'
	}

	// Clean up staged scalars
	scalar drop _aaivaft_ep_hls_l  _aaivaft_ep_hls_g  _aaivaft_ep_hls_e  ///
	            _aaivaft_ep_ci_l_l _aaivaft_ep_ci_u_l                   ///
	            _aaivaft_ep_ci_l_g _aaivaft_ep_ci_u_g                   ///
	            _aaivaft_ep_ci_l_e _aaivaft_ep_ci_u_e                   ///
	            _aaivaft_ep_nobs   _aaivaft_ep_has_inst                  ///
	            _aaivaft_ep_do_ls  _aaivaft_ep_do_rank _aaivaft_ep_level

	tempname b V

	// E(b) holds point estimates only; bounds go to e() scalars.
	// A zero V is posted because estimates table, outreg2 and esttab
	// require e(V); aaivaft gives confidence sets, not standard errors.
	// Pieces are added only when computed, since e(b) cannot hold missings.
	local bvals ""
	local bnames ""
	if `do_rank' {
		local bvals  "`hls_l', `hls_g', `hls_e'"
		local bnames "HLS_logistic HLS_gumbel HLS_exponential"
		if `has_fo' {
			local bvals  "`bvals', `hls_q'"
			local bnames "`bnames' HLS_empq"
		}
	}
	if `has_inst' & `do_ls' {
		if "`bvals'" != "" local bvals "`bvals', "
		local bvals  "`bvals'`ls_b'"
		local bnames "`bnames' LS_IV"
	}
	matrix `b' = (`bvals')
	matrix colnames `b' = `bnames'
	matrix rownames `b' = `treatvar'

	local k = colsof(`b')
	matrix `V' = J(`k', `k', 0)
	matrix rownames `V' = `bnames'
	matrix colnames `V' = `bnames'

	// Post to e()
	ereturn post `b' `V', obs(`nobs') depname(`depvar') esample(`touse')

	ereturn local  cmd         "aaivaft"
	ereturn local  depvar      "`depvar'"
	ereturn local  treatvar    "`treatvar'"
	ereturn local  title       "Rank-based IV AFT estimator"
	if `has_inst' & "`instrument'" != "" {
		ereturn local instrument "`instrument'"
	}
	ereturn scalar level         = `lvl'
	ereturn scalar N             = `nobs'
	ereturn scalar HLS_logistic  = `hls_l'
	ereturn scalar HLS_gumbel    = `hls_g'
	ereturn scalar HLS_exponential = `hls_e'

	// Rank confidence sets
	if `has_inst' & `do_rank' {
		ereturn scalar CIR_lo_logistic  = `cr_lo_l'
		ereturn scalar CIR_hi_logistic  = `cr_hi_l'
		ereturn scalar CIR_lo_gumbel    = `cr_lo_g'
		ereturn scalar CIR_hi_gumbel    = `cr_hi_g'
		ereturn scalar CIR_lo_exponential = `cr_lo_e'
		ereturn scalar CIR_hi_exponential = `cr_hi_e'
		if `has_fo' {
			ereturn local  fo          `"$AAIVAFT_ep_fo"'
			ereturn scalar HLS_empq    = `hls_q'
			ereturn scalar CIR_lo_empq = `cr_lo_q'
			ereturn scalar CIR_hi_empq = `cr_hi_q'
		}
	}

	if `has_inst' & `do_ls' {
		ereturn scalar CI_lo_logistic  = `ci_l_l'
		ereturn scalar CI_hi_logistic  = `ci_u_l'
		ereturn scalar CI_lo_gumbel    = `ci_l_g'
		ereturn scalar CI_hi_gumbel    = `ci_u_g'
		ereturn scalar CI_lo_exponential = `ci_l_e'
		ereturn scalar CI_hi_exponential = `ci_u_e'
		ereturn scalar beta_LS           = `ls_b'
	}

end




// =============================================================================
// Helper: projection P*mu via FWL
// =============================================================================
mata:
matrix _aaivaft_Pmu(matrix M1z, matrix Pz, matrix mu)
{
	return(M1z * (Pz * (M1z' * mu)))
}
end

// =============================================================================
// Helper: annihilator M*mu via FWL
// =============================================================================
mata:
matrix _aaivaft_Mmu(matrix x, matrix Qx, matrix mu)
{
	return(mu - x * (Qx * (x' * mu)))
}
end

// =============================================================================
// Helper: quadratic inversion for confidence set bounds
// When a < 0, the two roots are swapped — sort them
// =============================================================================
mata:
matrix _aaivaft_ci(real scalar Tc,
                 real scalar k2, real scalar n,  real scalar k,
                 real scalar YPY, real scalar yPY, real scalar yPy,
                 real scalar YMY, real scalar yMY, real scalar yMy)
{
	G   = (Tc * k2) / (n - k)
	a   =    YPY - G * YMY
	b   = -2*(yPY - G * yMY)
	c   =    yPy - G * yMy
	det = b^2 - 4*a*c
	if (det < 0) return(. \ .)

	r1 = (-b - sqrt(det)) / (2*a)
	r2 = (-b + sqrt(det)) / (2*a)

	if (r1 <= r2) return(r1 \ r2)
	else          return(r2 \ r1)
}
end


// =============================================================================
// Helper: rank confidence set on the grid
// Returns (lower, upper, edge flag, gap flag). Keeps grid values with G < gc.
// =============================================================================
mata:
real rowvector _aaivaft_rankcs(real colvector grid, real colvector G, real scalar gc)
{
	real colvector acc
	real scalar    m, lo, hi, edge, gap

	acc = selectindex(G :< gc)
	m   = rows(acc)
	if (m == 0) return((., ., 0, 0))

	lo   = grid[acc[1]]
	hi   = grid[acc[m]]
	edge = (acc[1] == 1 | acc[m] == rows(grid))
	gap  = (m < acc[m] - acc[1] + 1)
	return((lo, hi, edge, gap))
}
end


// =============================================================================
// Helper: empirical quantile score Fo^-1(i/(n+1)), i = 1..n, from the
// draws d, with linear interpolation between order statistics
// =============================================================================
mata:
real colvector _aaivaft_foquant(real colvector d, real scalar n)
{
	real colvector s, u, h, lo, w
	real scalar    B

	s  = sort(d, 1)
	B  = rows(s)
	u  = (1::n) :/ (n + 1)
	h  = u :* (B - 1) :+ 1
	lo = rowmin((floor(h), J(n, 1, B - 1)))
	w  = h - lo
	return(s[lo] :* (1 :- w) + s[lo :+ 1] :* w)
}
end


// =============================================================================
// Mata function 1a: Rank IV — with instrument (original)
// =============================================================================
mata:
void _aaivaft_rankinf(
	string scalar sY,
	string scalar sT,
	string scalar sZ,
	string scalar sX,
	string scalar sgmin,
	string scalar sgmax,
	string scalar sgstep,
	string scalar sreps,
	string scalar slevel,
	real scalar   hasfo,
	real scalar   zres
)
{
	y  = st_matrix(sY)
	Y  = st_matrix(sT)
	z  = st_matrix(sZ)
	x1 = st_matrix(sX)
	n  = rows(y)

	gmin  = st_numscalar(sgmin)
	gmax  = st_numscalar(sgmax)
	gstep = st_numscalar(sgstep)
	grid  = range(gmin, gmax, gstep)
	nb    = rows(grid)

	x1d = x1 :- mean(x1)
	// FWL form, avoids building the n x n matrix I(n) - x1d*Q*x1d'
	Q1d = invsym(x1d' * x1d)
	M1y = y - x1d * (Q1d * (x1d' * y))
	M1Y = Y - x1d * (Q1d * (x1d' * Y))

	// By default use the instrument residualized on the
	// covariates, M1*z; with zcenter only center it
	if (zres) zd = (z :- mean(z)) - x1d * (Q1d * (x1d' * z))
	else      zd = z :- mean(z)
	zdzd = (zd' * zd)[1,1]

	// Deterministic tie-breaking: tiny jitter proportional to obs index
	// Needed when treatment is discrete (e.g. integer severity scores)
	// Magnitude 1e-12 is negligible relative to any real variation in E
	jitter = (1::n) :* (1e-12)

	// Empirical quantile score at ranks 1..n
	if (hasfo) {
		pfo = findexternal("_aaivaft_fodraws")
		s_q = _aaivaft_foquant(*pfo, n)
	}

	Rx = J(nb, 4 + hasfo, .)

	for (i = 1; i <= nb; i++) {
		E  = M1y - M1Y * grid[i] + jitter
		R  = mm_ranks(E, 1, 0)
		Rp = R :/ (n + 1)
		if (hasfo) {
			c_q = s_q[round(R)]
			Rx[i, 5] = ((zd' * c_q)[1,1])^2 / zdzd
		}

		c_l = ln(Rp :/ (1 :- Rp))
		c_g = -ln(-ln(Rp))
		c_e = -ln(1 :- Rp)

		// C_l, c_g, c_e are in data order, so they pair with zd directly
		zdcl = (zd' * c_l)[1,1]
		zdcg = (zd' * c_g)[1,1]
		zdce = (zd' * c_e)[1,1]

		Rx[i, 1] = grid[i]
		Rx[i, 2] = zdcl^2 / zdzd
		Rx[i, 3] = zdcg^2 / zdzd
		Rx[i, 4] = zdce^2 / zdzd
	}

	// Use direct minimum search — robust to ties (mm_ranks+selectindex fails on ties)
	idx_l  = selectindex(Rx[.,2] :== min(Rx[.,2]))
	idx_g  = selectindex(Rx[.,3] :== min(Rx[.,3]))
	idx_e  = selectindex(Rx[.,4] :== min(Rx[.,4]))

	beta_l = mean(Rx[idx_l, 1])
	beta_g = mean(Rx[idx_g, 1])
	beta_e = mean(Rx[idx_e, 1])

	st_numscalar("r(beta_logistic)",  beta_l)
	st_numscalar("r(beta_gumbel)",    beta_g)
	st_numscalar("r(beta_exponential)", beta_e)
	if (hasfo) {
		idx_q  = selectindex(Rx[.,5] :== min(Rx[.,5]))
		st_numscalar("r(beta_empq)", mean(Rx[idx_q, 1]))
	}

	// -------------------------------------------------------------------------
	// Rank confidence sets by simulation
	// Under H0: beta = b0 the ranks of E are a random permutation, so the null
	// distribution of G does not depend on b0. Pairing a permuted zd with the
	// scores at ranks 1..n gives the same distribution as drawing U ~ U(0,1)^n,
	// ranking it and applying the score.
	// -------------------------------------------------------------------------
	reps  = st_numscalar(sreps)
	level = st_numscalar(slevel)
	idx   = floor(reps * level / 100)

	u0  = (1::n) :/ (n + 1)
	s_l = ln(u0 :/ (1 :- u0))
	s_g = -ln(-ln(u0))
	s_e = -ln(1 :- u0)

	// Permute a sorted copy so the draws do not depend on data order
	zsrt = sort(zd, 1)
	sim = J(reps, 3 + hasfo, .)
	for (r = 1; r <= reps; r++) {
		zp = zsrt[order(runiform(n, 1), 1)]
		sim[r, 1] = (zp' * s_l)^2 / zdzd
		sim[r, 2] = (zp' * s_g)^2 / zdzd
		sim[r, 3] = (zp' * s_e)^2 / zdzd
		if (hasfo) sim[r, 4] = (zp' * s_q)^2 / zdzd
	}

	gc_l = sort(sim[., 1], 1)[idx]
	gc_g = sort(sim[., 2], 1)[idx]
	gc_e = sort(sim[., 3], 1)[idx]

	cs_l = _aaivaft_rankcs(Rx[., 1], Rx[., 2], gc_l)
	cs_g = _aaivaft_rankcs(Rx[., 1], Rx[., 3], gc_g)
	cs_e = _aaivaft_rankcs(Rx[., 1], Rx[., 4], gc_e)

	st_numscalar("r(cir_l_logistic)",     cs_l[1])
	st_numscalar("r(cir_u_logistic)",     cs_l[2])
	st_numscalar("r(cir_edge_logistic)",  cs_l[3])
	st_numscalar("r(cir_gap_logistic)",   cs_l[4])
	st_numscalar("r(gcrit_logistic)",     gc_l)
	st_numscalar("r(cir_l_gumbel)",       cs_g[1])
	st_numscalar("r(cir_u_gumbel)",       cs_g[2])
	st_numscalar("r(cir_edge_gumbel)",    cs_g[3])
	st_numscalar("r(cir_gap_gumbel)",     cs_g[4])
	st_numscalar("r(gcrit_gumbel)",       gc_g)
	st_numscalar("r(cir_l_exponential)",    cs_e[1])
	st_numscalar("r(cir_u_exponential)",    cs_e[2])
	st_numscalar("r(cir_edge_exponential)", cs_e[3])
	st_numscalar("r(cir_gap_exponential)",  cs_e[4])
	st_numscalar("r(gcrit_exponential)",    gc_e)

	if (hasfo) {
		gc_q = sort(sim[., 4], 1)[idx]
		cs_q = _aaivaft_rankcs(Rx[., 1], Rx[., 5], gc_q)
		st_numscalar("r(cir_l_empq)",    cs_q[1])
		st_numscalar("r(cir_u_empq)",    cs_q[2])
		st_numscalar("r(cir_edge_empq)", cs_q[3])
		st_numscalar("r(cir_gap_empq)",  cs_q[4])
		st_numscalar("r(gcrit_empq)",    gc_q)
	}
}
end


// =============================================================================
// Mata function 1b: Rank — WITHOUT instrument (rank-only benchmark)
// G_RO(b0) = (s~' c)^2 / (s~' s~), with s~ = M1*s the treatment
// residualized on the covariates, as in the paper
// =============================================================================
mata:
void _aaivaft_rankinf_noinst(
	string scalar sY,
	string scalar sT,
	string scalar sX,
	string scalar sgmin,
	string scalar sgmax,
	string scalar sgstep
)
{
	y  = st_matrix(sY)
	Y  = st_matrix(sT)
	x1 = st_matrix(sX)
	n  = rows(y)

	gmin  = st_numscalar(sgmin)
	gmax  = st_numscalar(sgmax)
	gstep = st_numscalar(sgstep)
	grid  = range(gmin, gmax, gstep)
	nb    = rows(grid)

	// Residualise y and Y on covariates
	x1d = x1 :- mean(x1)
	// FWL form, avoids building the n x n matrix I(n) - x1d*Q*x1d'
	Q1d = invsym(x1d' * x1d)
	M1y = y - x1d * (Q1d * (x1d' * y))
	M1Y = Y - x1d * (Q1d * (x1d' * Y))

	// Residualized treatment (fixed across the beta grid)
	RYdd = M1Y :- mean(M1Y)
	RYss = (RYdd' * RYdd)[1,1]

	// Deterministic tie-breaking (same rationale as _aaivaft_rankinf)
	jitter = (1::n) :* (1e-12)

	Rx = J(nb, 4, .)

	for (i = 1; i <= nb; i++) {
		E  = M1y - M1Y * grid[i] + jitter
		R  = mm_ranks(E, 1, 0)
		Rp = R :/ (n + 1)

		c_l = ln(Rp :/ (1 :- Rp))
		c_g = -ln(-ln(Rp))
		c_e = -ln(1 :- Rp)

		// scores are in data order, so they pair with RYdd directly
		RYcl = (RYdd' * c_l)[1,1]
		RYcg = (RYdd' * c_g)[1,1]
		RYce = (RYdd' * c_e)[1,1]

		Rx[i, 1] = grid[i]
		Rx[i, 2] = RYcl^2 / RYss
		Rx[i, 3] = RYcg^2 / RYss
		Rx[i, 4] = RYce^2 / RYss
	}

	// Use direct minimum search — robust to ties
	idx_l  = selectindex(Rx[.,2] :== min(Rx[.,2]))
	idx_g  = selectindex(Rx[.,3] :== min(Rx[.,3]))
	idx_e  = selectindex(Rx[.,4] :== min(Rx[.,4]))

	beta_l = mean(Rx[idx_l, 1])
	beta_g = mean(Rx[idx_g, 1])
	beta_e = mean(Rx[idx_e, 1])

	st_numscalar("r(beta_logistic)",  beta_l)
	st_numscalar("r(beta_gumbel)",    beta_g)
	st_numscalar("r(beta_exponential)", beta_e)
}
end


// =============================================================================
// Mata function 2: Least-Squares IV  [FAST] — requires instrument
// =============================================================================
mata:
void _aaivaft_leastsquares(
	string scalar sY,
	string scalar sT,
	string scalar sZ,
	string scalar sX,
	string scalar sreps,
	string scalar slevel
)
{
	y  = st_matrix(sY)
	Y  = st_matrix(sT)
	z  = st_matrix(sZ)
	x1 = st_matrix(sX)

	// Add a constant; k1 is the rank, so a full set of dummies
	// plus the constant does not inflate the degrees of freedom
	x1 = J(rows(x1), 1, 1), x1

	n     = rows(y)
	k1    = rank(x1)
	k2    = cols(z)
	k     = k1 + k2
	reps  = st_numscalar(sreps)
	level = st_numscalar(slevel)
	alpha = 1 - level / 100
	idx   = floor(reps * (1 - alpha))

	x   = x1, z
	Qx  = invsym(x' * x)
	Qx1 = invsym(x1' * x1)

	M1z = z - x1 * (Qx1 * (x1' * z))
	Pz  = invsym(M1z' * M1z)

	// Simulate on rows sorted by (z, covariates), so the draws do not
	// depend on data order; Qx and Pz do not change under row reordering
	ps   = order((z, x1), (1..(cols(z) + cols(x1))))
	M1zs = M1z[ps, .]
	xs   = x[ps, .]

	U    = runiform(n, reps)
	mu_l = ln(U  :/ (1 :- U))
	mu_g = -ln(-ln(U))
	mu_e = -ln(1 :- U)

	Pmu_l = _aaivaft_Pmu(M1zs, Pz, mu_l)
	Mmu_l = _aaivaft_Mmu(xs, Qx, mu_l)
	GAR_l = (colsum(Pmu_l :* mu_l) :/ colsum(Mmu_l :* mu_l)) :* ((n-k)/k2)

	Pmu_g = _aaivaft_Pmu(M1zs, Pz, mu_g)
	Mmu_g = _aaivaft_Mmu(xs, Qx, mu_g)
	GAR_g = (colsum(Pmu_g :* mu_g) :/ colsum(Mmu_g :* mu_g)) :* ((n-k)/k2)

	Pmu_e = _aaivaft_Pmu(M1zs, Pz, mu_e)
	Mmu_e = _aaivaft_Mmu(xs, Qx, mu_e)
	GAR_e = (colsum(Pmu_e :* mu_e) :/ colsum(Mmu_e :* mu_e)) :* ((n-k)/k2)

	Tc_l = sort(GAR_l', 1)[idx, 1]
	Tc_g = sort(GAR_g', 1)[idx, 1]
	Tc_e = sort(GAR_e', 1)[idx, 1]

	M1zY = M1z' * Y
	M1zy = M1z' * y

	YPY = (M1zY' * Pz * M1zY)[1,1]
	yPY = (M1zy' * Pz * M1zY)[1,1]
	yPy = (M1zy' * Pz * M1zy)[1,1]

	xY  = x' * Y
	xy  = x' * y
	YMY = (Y'*Y  - xY'*Qx*xY)[1,1]
	yMY = (y'*Y  - xy'*Qx*xY)[1,1]
	yMy = (y'*y  - xy'*Qx*xy)[1,1]

	// LS IV point estimate, the value of beta where GAR is zero
	beta_ls = yPY / YPY

	ci_l = _aaivaft_ci(Tc_l, k2, n, k, YPY, yPY, yPy, YMY, yMY, yMy)
	ci_g = _aaivaft_ci(Tc_g, k2, n, k, YPY, yPY, yPy, YMY, yMY, yMy)
	ci_e = _aaivaft_ci(Tc_e, k2, n, k, YPY, yPY, yPy, YMY, yMY, yMy)

	st_numscalar("r(ci_l_logistic)",   ci_l[1])
	st_numscalar("r(ci_u_logistic)",   ci_l[2])
	st_numscalar("r(ci_l_gumbel)",     ci_g[1])
	st_numscalar("r(ci_u_gumbel)",     ci_g[2])
	st_numscalar("r(ci_l_exponential)",  ci_e[1])
	st_numscalar("r(ci_u_exponential)",  ci_e[2])
	st_numscalar("r(beta_ls)",           beta_ls)
}
end

