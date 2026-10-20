*! version 1.0.0  19oct2026
*! aaivaft_export: writes aaivaft results to a CSV file
*! Companion to aaivaft.ado; called by aaivaft's export() option or run by itself
*! after aaivaft. Syntax: aaivaft_export using filename[.csv] [, replace]

// =============================================================================
// aaivaft_export: export aaivaft results to a CSV file
// Syntax: aaivaft_export using filename[.csv] [, replace]
// Reads from e() after aaivaft. One row per score with the HLS point, the rank
// confidence set, the LS point and the LS confidence set. Empty cells mark
// results that were not computed. A name without .csv gets the .csv suffix.
// =============================================================================
program define aaivaft_export
    version 14.0

    syntax using/ [, replace]

    if "`e(cmd)'" != "aaivaft" {
        di as error "Last estimation command was not aaivaft."
        exit 301
    }

    // Force the .csv suffix
    mata: st_local("suf", pathsuffix(st_local("using")))
    if "`suf'" == "" {
        local using "`using'.csv"
    }
    else if lower("`suf'") != ".csv" {
        mata: st_local("using", pathrmsuffix(st_local("using")) + ".csv")
        di as txt " Note: aaivaft_export writes CSV; file saved as " as res `""`using'""'
    }

    capture confirm file `"`using'"'
    local exists = (_rc == 0)
    if `exists' & "`replace'" == "" {
        di as error "File `using' already exists. Use replace option."
        exit 602
    }
    // Use replace only when the file exists, to avoid Stata's
    // "file not found" note for a new file
    local fopt ""
    if `exists' local fopt "replace"

    local depvar   "`e(depvar)'"
    local treatvar "`e(treatvar)'"
    local inst     "`e(instrument)'"
    local nobs     = e(N)
    local lvl      = e(level)

    tempname fh
    file open `fh' using `"`using'"', write text `fopt'
    file write `fh' "score,hls,rank_lo,rank_hi,ls_point,ls_lo,ls_hi,"
    file write `fh' "depvar,treatvar,instrument,N,level" _n

    // Add the empirical quantile row when aaivaft was run with fo()
    local scores "logistic gumbel exponential"
    if !missing(e(HLS_empq)) local scores "`scores' empq"

    foreach s of local scores {
        local row "`s'"
        foreach q in HLS_`s' CIR_lo_`s' CIR_hi_`s' beta_LS CI_lo_`s' CI_hi_`s' {
            local v = e(`q')
            if missing(`v') {
                local row "`row',"
            }
            else {
                local sv = strtrim(string(`v', "%12.6f"))
                local row "`row',`sv'"
            }
        }
        file write `fh' `"`row',`depvar',`treatvar',`inst',`nobs',`lvl'"' _n
    }
    file close `fh'

    di as txt " Results exported to: " as res `""`using'""'

end
