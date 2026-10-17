*! version 0.1.3  17oct2026  powerape for Stata -- power analysis for average
*! partial effects (APE) and average interaction effects (AIE). Port of the
*! powerape R package 1.11.3 (Jesper Wulff). Estimation inside the
*! simulation loop uses Stata's own probit/logit/cfprobit (engine(stata), the
*! default) with closed-form delta-method effects validated against margins
*! and ginteff, or margins/ginteff themselves (effects(margins)).
program powerape, rclass
    version 16
    gettoken sub 0 : 0, parse(" ,")
    local sub = lower(`"`sub'"')
    local known dgp setape setaie describe list drop trueape trueaie power ///
        curve n mde robust alpha powerat frontier statement engine version _parse
    if `"`sub'"' == "" {
        di as err "powerape: a subcommand is required (" `"`known'"' ")"
        exit 198
    }
    if !`: list sub in known' {
        di as err `"powerape: unknown subcommand "`sub'""'
        exit 198
    }
    if "`sub'" == "n" local sub nreq
    powerape_`sub' `0'
    return add
end

program powerape__parse, sclass
    pape_parsevars `0'
    sreturn list
end

* ---------------------------------------------------------------------------
* helpers
* ---------------------------------------------------------------------------

* parse one variable spec: "name[, binary|normal p(#) mean(#) sd(#) icc(#)]"
program pape_parseone, sclass
    sreturn clear
    syntax name(name=vname id="variable name") [, Binary Normal p(string) ///
        mean(string) sd(string) icc(string)]
    local type ""
    if "`binary'" != "" local type binary
    if "`normal'" != "" {
        if "`type'" != "" {
            di as err "variable `vname': choose binary or normal, not both"
            exit 198
        }
        local type normal
    }
    if "`type'" == "" {
        if "`p'" != "" local type binary
        else local type normal
    }
    foreach o in p mean sd icc {
        if "``o''" == "" local `o' .
        else {
            capture confirm number ``o''
            if _rc {
                di as err "variable `vname': `o'(``o'') is not a number"
                exit 198
            }
        }
    }
    sreturn local name `vname'
    sreturn local type `type'
    sreturn local p `p'
    sreturn local mean `mean'
    sreturn local sd `sd'
    sreturn local icc `icc'
end

* parse a list of specs "(spec) (spec)" or a single spec into parallel lists
program pape_parsevars, sclass
    args specs
    sreturn clear
    local names ""
    local types ""
    local ps ""
    local means ""
    local sds ""
    local iccs ""
    local rest = strtrim(`"`specs'"')
    if substr(`"`rest'"', 1, 1) != "(" & `"`rest'"' != "" {
        local rest `"(`rest')"'
    }
    while `"`rest'"' != "" {
        gettoken spec rest : rest, match(par)
        if `"`par'"' != "(" {
            di as err `"could not parse variable specification: `spec'"'
            exit 198
        }
        pape_parseone `spec'
        local names `names' `s(name)'
        local types `types' `s(type)'
        local ps `ps' `s(p)'
        local means `means' `s(mean)'
        local sds `sds' `s(sd)'
        local iccs `iccs' `s(icc)'
        local rest = strtrim(`"`rest'"')
    }
    sreturn local names `names'
    sreturn local types `types'
    sreturn local ps `ps'
    sreturn local means `means'
    sreturn local sds `sds'
    sreturn local iccs `iccs'
    local k : word count `names'
    sreturn local k `k'
end

* correlation option: "" | scalar | matrix name
program pape_parsecorr, sclass
    args corr
    sreturn clear
    local corr = strtrim(`"`corr'"')
    if `"`corr'"' == "" {
        sreturn local kind 0
        sreturn local scalar .
        sreturn local matname ""
        exit
    }
    capture confirm number `corr'
    if !_rc {
        sreturn local kind 1
        sreturn local scalar `corr'
        sreturn local matname ""
        exit
    }
    capture confirm matrix `corr'
    if _rc {
        di as err "correlation() must be a number or the name of a matrix"
        exit 198
    }
    sreturn local kind 2
    sreturn local scalar .
    sreturn local matname `corr'
end

* set engine defaults for this call
program pape_setengine
    syntax [, ENGine(string) EFFects(string)]
    * remember the session default (set by -powerape engine-) so that
    * pape_runmata restores it afterwards: a per-command engine()/effects()
    * applies to that command only (2026-10-16; before, it leaked into the
    * session default)
    mata: st_global("PAPE_ENG_PREV", pape_engine())
    mata: st_global("PAPE_EFF_PREV", pape_effects())
    if "`engine'" == "" local engine $PAPE_ENG_PREV
    if "`effects'" == "" local effects $PAPE_EFF_PREV
    local engine = lower("`engine'")
    local effects = lower("`effects'")
    if "`engine'" == "closed" | "`engine'" == "mata" local engine = cond("`engine'" == "closed", "stata", "mata")
    mata: pape_set_engine("`engine'", "`effects'")
end

* run a Mata command with frame safety
program pape_runmata
    args cmd restore
    local userframe = c(frame)
    * the Stata engine fits probit/logit/cfprobit, which would replace the
    * user's estimation results: hold them (restored even on error)
    tempname eh
    capture _estimates hold `eh', restore
    capture noisily mata: `cmd'
    local rc = _rc
    capture _estimates unhold `eh'
    capture mata: pape_set_quiet(0)
    * the simulation commands pass -restore-: put back the session engine
    * default that pape_setengine saved (not done for the nested -powerape dgp-
    * calls a sweep issues while its own engine choice must stay in force)
    if "`restore'" != "" & `"$PAPE_ENG_PREV"' != "" {
        capture mata: pape_set_engine("$PAPE_ENG_PREV", "$PAPE_EFF_PREV")
        macro drop PAPE_ENG_PREV PAPE_EFF_PREV
    }
    capture frame change `userframe'
    capture frame drop _pape_sim
    if `rc' exit `rc'
end


* ---------------------------------------------------------------------------
* powerape dgp
* ---------------------------------------------------------------------------
program powerape_dgp, rclass
    syntax name(name=dname id="DGP name") [if] [in] [, EMPirical FROMfit PANel IV ///
        Model(string) Focal(string) MODerator(string) COVariates(string) ///
        INSTruments(string) CORRelation(string) BASEline(string) SIGnal(real 0) ///
        GAMma(numlist) NINT(real -1) SEEDint(integer 20260812) ///
        NPERiods(numlist integer) PPERiods(numlist) RETention(string) ///
        RHO(real 0) CREshare(real 0) ENDOGeneity(string) IVStrength(string) ///
        REPLACE noDISplay]
    local spec_opts : copy local 0
    gettoken junk spec_opts : spec_opts, parse(",")
    local spec_opts = subinstr(`"`spec_opts'"', ",", "", 1)
    local spec_opts = subinstr(`"`spec_opts'"', " replace", "", .)
    local spec_opts = subinstr(`"`spec_opts'"', " nodisplay", "", .)
    if "`model'" == "" local model probit
    local model = lower("`model'")
    local disp = ("`display'" == "")
    local nroutes = ("`empirical'" != "") + ("`fromfit'" != "") + ("`panel'" != "") + ("`iv'" != "")
    if `nroutes' > 1 {
        di as err "choose at most one of empirical, fromfit, panel, iv"
        exit 198
    }
    if "`replace'" == "" {
        mata: st_local("exists", strofreal(pape_store_has("DGP", "`dname'")))
        if `exists' {
            di as err "DGP `dname' already exists; use replace"
            exit 110
        }
    }
    if "`baseline'" != "" {
        capture confirm number `baseline'
        if _rc {
            di as err "baseline() must be a number"
            exit 198
        }
    }
    pape_parsecorr `"`correlation'"'
    local ckind `s(kind)'
    local cscalar `s(scalar)'
    local cmat `s(matname)'

    * ---- fromfit -----------------------------------------------------------
    if "`fromfit'" != "" {
        if "`focal'" == "" {
            di as err "fromfit needs focal(varname)"
            exit 198
        }
        if "`baseline'" == "" local baseline .
        tempvar touse
        quietly gen byte `touse' = e(sample)
        if "`if'`in'" != "" {
            quietly replace `touse' = 0 if !(`touse')
            marksample touse2
            quietly replace `touse' = `touse' & `touse2'
        }
        mata: pape_cmd_dgp_fromfit("`dname'", "`focal'", `baseline', "`touse'", `"`spec_opts'"', `disp')
        pape_return_dgp
        return add
        exit
    }

    if "`baseline'" == "" {
        di as err "baseline() is required"
        exit 198
    }

    * ---- empirical ---------------------------------------------------------
    if "`empirical'" != "" {
        marksample touse, novarlist
        if "`focal'" == "" {
            di as err "empirical needs focal(varname) or focal(name, binary|normal ...)"
            exit 198
        }
        local focal = strtrim(`"`focal'"')
        if strpos(`"`focal'"', ",") == 0 & `: word count `focal'' == 1 {
            capture confirm variable `focal'
            if _rc {
                di as err "focal(`focal') is not a variable in the data (an independent focal needs a spec with a type)"
                exit 198
            }
            local joint 1
            local varlist `focal' `covariates'
            local fname ""
            local ftype ""
            local fp .
            local fmean .
            local fsd .
        }
        else {
            local joint 0
            pape_parsevars `"`focal'"'
            if `s(k)' != 1 {
                di as err "focal() must describe one variable"
                exit 198
            }
            local fname `s(names)'
            local ftype `s(types)'
            local fp `s(ps)'
            local fmean `s(means)'
            local fsd `s(sds)'
            local varlist `covariates'
        }
        if "`covariates'" != "" {
            capture confirm variable `covariates'
            if _rc {
                di as err "covariates() must list variables in the data for the empirical route"
                exit 198
            }
        }
        local varlist : list uniq varlist
        markout `touse' `varlist'
        local gstr "`gamma'"
        mata: pape_cmd_dgp_emp("`dname'", "`model'", "`varlist'", "`touse'", `joint', ///
            "`fname'", "`ftype'", `fp', `fmean', `fsd', `baseline', `signal', "`gstr'", ///
            `seedint', `"`spec_opts'"', `disp')
        pape_return_dgp
        return add
        exit
    }

    * ---- parametric / panel / iv -------------------------------------------
    if "`focal'" == "" {
        di as err "focal() is required"
        exit 198
    }
    pape_parsevars `"`focal'"'
    if `s(k)' != 1 {
        di as err "focal() must describe one variable"
        exit 198
    }
    local vnames `s(names)'
    local vtypes `s(types)'
    local vps `s(ps)'
    local vmeans `s(means)'
    local vsds `s(sds)'
    local viccs `s(iccs)'
    local has_m 0
    if `"`moderator'"' != "" {
        pape_parsevars `"`moderator'"'
        if `s(k)' != 1 {
            di as err "moderator() must describe one variable"
            exit 198
        }
        local has_m 1
        local vnames `vnames' `s(names)'
        local vtypes `vtypes' `s(types)'
        local vps `vps' `s(ps)'
        local vmeans `vmeans' `s(means)'
        local vsds `vsds' `s(sds)'
        local viccs `viccs' `s(iccs)'
    }
    if `"`covariates'"' != "" {
        pape_parsevars `"`covariates'"'
        local vnames `vnames' `s(names)'
        local vtypes `vtypes' `s(types)'
        local vps `vps' `s(ps)'
        local vmeans `vmeans' `s(means)'
        local vsds `vsds' `s(sds)'
        local viccs `viccs' `s(iccs)'
    }
    if "`panel'" != "" {
        if "`nperiods'" == "" {
            di as err "panel needs nperiods()"
            exit 198
        }
        if "`retention'" == "" local retention .
        if `nint' < 0 local nint 1e5
        if `seedint' == 20260812 local seedint 20260814
        mata: pape_cmd_dgp_panel("`dname'", "`model'", "`vnames'", "`vtypes'", "`vps'", ///
            "`vmeans'", "`vsds'", "`viccs'", `has_m', "`nperiods'", "`pperiods'", `retention', ///
            `rho', `creshare', `ckind', `cscalar', "`cmat'", `baseline', `signal', `nint', ///
            `seedint', `"`spec_opts'"', `disp')
        pape_return_dgp
        return add
        exit
    }
    if "`iv'" != "" {
        if `"`instruments'"' == "" {
            di as err "iv needs instruments()"
            exit 198
        }
        if "`endogeneity'" == "" | "`ivstrength'" == "" {
            di as err "iv needs endogeneity() and ivstrength()"
            exit 198
        }
        pape_parsevars `"`instruments'"'
        local inames `s(names)'
        local itypes `s(types)'
        local ips `s(ps)'
        local imeans `s(means)'
        local isds `s(sds)'
        if `nint' < 0 local nint 1e5
        if `seedint' == 20260812 local seedint 20260814
        mata: pape_cmd_dgp_iv("`dname'", "`vnames'", "`vtypes'", "`vps'", "`vmeans'", "`vsds'", ///
            "`viccs'", `has_m', "`inames'", "`itypes'", "`ips'", "`imeans'", "`isds'", ///
            `endogeneity', `ivstrength', `ckind', `cscalar', "`cmat'", `baseline', `signal', ///
            `nint', `seedint', `"`spec_opts'"', `disp')
        pape_return_dgp
        return add
        exit
    }
    if `nint' < 0 local nint 4e5
    local gstr "`gamma'"
    mata: pape_cmd_dgp_par("`dname'", "`model'", "`vnames'", "`vtypes'", "`vps'", "`vmeans'", ///
        "`vsds'", "`viccs'", `has_m', `ckind', `cscalar', "`cmat'", `baseline', `signal', ///
        "`gstr'", `nint', `seedint', `"`spec_opts'"', `disp')
    pape_return_dgp
    return add
end

program pape_return_dgp, rclass
    foreach v in model route estimand focal {
        return local `v' `"${pape_r_`v'}"'
    }
    foreach v in baseline signal beta0 beta_focal beta_mod beta_int target p1 k {
        capture return scalar `v' = scalar(pape_r_`v')
    }
    capture return matrix gamma = pape_r_gamma
    pape_clear_r
end

program pape_clear_r
    capture macro drop pape_r_*
    foreach v in baseline signal beta0 beta_focal beta_mod beta_int target p1 k power mcse n nsim alpha conf sesoi separated min_cell n_failed goal confirmed nsim_confirm mde is_margin value worst nmax nmax_power nmax_mcse nmax_scenario beta wcer power_default wcer_default range_lo range_hi capped cap alpha_uncapped wcer_uncapped size size_mcse size_truth {
        capture scalar drop pape_r_`v'
    }
end

* ---------------------------------------------------------------------------
* setape / setaie / describe / trueape / trueaie / list / drop / version / engine
* ---------------------------------------------------------------------------
program powerape_setape, rclass
    syntax name(name=dname id="DGP name"), Target(real) [noDISplay]
    local disp = ("`display'" == "")
    mata: pape_cmd_setape("`dname'", `target', `disp')
    pape_return_dgp
    return add
end

program powerape_setaie, rclass
    syntax name(name=dname id="DGP name"), Target(real) MAINFocal(real) MAINModerator(real) [noDISplay]
    local disp = ("`display'" == "")
    mata: pape_cmd_setaie("`dname'", `target', `mainfocal', `mainmoderator', `disp')
    pape_return_dgp
    return add
end

program powerape_describe, rclass
    syntax name(name=dname id="DGP name")
    mata: pape_cmd_describe("`dname'")
    pape_return_dgp
    return add
end

program powerape_trueape, rclass
    syntax name(name=dname id="DGP name")
    mata: pape_cmd_true("`dname'", "ape")
    return scalar ape = scalar(pape_r_value)
    pape_clear_r
end

program powerape_trueaie, rclass
    syntax name(name=dname id="DGP name")
    mata: pape_cmd_true("`dname'", "aie")
    return scalar aie = scalar(pape_r_value)
    pape_clear_r
end

program powerape_list, rclass
    mata: pape_cmd_list()
    return local dgps `"${pape_r_dgps}"'
    return local results `"${pape_r_results}"'
    pape_clear_r
end

program powerape_drop
    syntax name(name=dname id="name")
    mata: pape_cmd_drop("`dname'")
end

program powerape_version, rclass
    * the version string of this ado-file; must equal the first word of
    * pape_version() in lpowerape.mlib (powerape_search.mata)
    local adov "0.1.3"
    capture mata: st_local("v", pape_version())
    if _rc {
        local rc = _rc
        di as err "powerape: the Mata library lpowerape.mlib is not available (pape_version() not found, rc `rc')."
        di as err "  After installing or rebuilding the library, run:  mata: mata mlib index"
        di as err "  The library shipped with the package was compiled by StataNow 19.5; on Stata 16-18"
        di as err "  rebuild it from the .mata sources (net get powerape, then do build.do)."
        exit `rc'
    }
    local mv : word 1 of `v'
    quietly findfile powerape.ado
    local adofile `"`r(fn)'"'
    capture quietly findfile lpowerape.mlib
    if _rc local mlibfile "(not found on the adopath; the functions are loaded in memory)"
    else local mlibfile `"`r(fn)'"'
    di as txt "powerape for Stata, version `v'"
    di as txt "  ado-file:      `adofile'"
    di as txt "  Mata library:  `mlibfile'"
    if "`mv'" != "`adov'" {
        di as err "warning: powerape.ado is version `adov' but lpowerape.mlib is version `mv'."
        di as err "  Rebuild the library with build.do (or reinstall the package) and run mata: mata mlib index."
    }
    else di as txt "  the ado-file and the Mata library agree (version `adov')"
    return scalar match = ("`mv'" == "`adov'")
    return local mlib `"`mlibfile'"'
    return local ado `"`adofile'"'
    return local mlib_version "`mv'"
    return local ado_version "`adov'"
    return local version "`v'"
end

program powerape_engine, rclass
    syntax [, ENGine(string) EFFects(string)]
    if "`engine'" == "" & "`effects'" == "" {
        mata: st_local("e", pape_engine()); mata: st_local("f", pape_effects())
        di as txt "engine(`e') effects(`f')"
        return local engine "`e'"
        return local effects "`f'"
        exit
    }
    mata: st_local("e", pape_engine()); mata: st_local("f", pape_effects())
    if "`engine'" == "" local engine `e'
    if "`effects'" == "" local effects `f'
    mata: pape_set_engine("`engine'", "`effects'")
    di as txt "engine(`engine') effects(`effects') set as the session default"
    return local engine "`engine'"
    return local effects "`effects'"
end


* ---------------------------------------------------------------------------
* common option parsing for the simulation commands
* ---------------------------------------------------------------------------
program pape_simopts, sclass
    syntax, [CLaim(string) SESoi(string) ALPha(string) ALPHARule(string) CONF(string) ///
        NSIM(real -1) SEED(integer -1) SE(string) SEParation(string) Name(string) ///
        DEFAULTclaim(string) DEFAULTnsim(real 1000) noDISplay noDRAWs ENGine(string) EFFects(string)]
    sreturn clear
    if "`claim'" == "" local claim `defaultclaim'
    local claim = lower("`claim'")
    if !inlist("`claim'", "minimum", "detect", "equivalence") {
        di as err "claim() must be minimum, detect, or equivalence"
        exit 198
    }
    if "`sesoi'" == "" local sesoi .
    if "`alpha'" == "" local alpha .
    if "`conf'" == "" local conf .
    foreach o in sesoi alpha conf {
        capture confirm number ``o''
        if _rc & "``o''" != "." {
            di as err "`o'() must be a number"
            exit 198
        }
    }
    if "`alpharule'" != "" {
        if "`alpha'" != "." {
            di as err "specify alpha() or alpharule(), not both"
            exit 198
        }
        if strpos(`"`alpharule'"', "@n") == 0 {
            di as err "alpharule() must be a Stata expression in @n, e.g. alpharule(0.05*sqrt(1000/@n))"
            exit 198
        }
    }
    if "`alpha'" == "." & "`alpharule'" == "" local alpha 0.05
    if `nsim' < 0 local nsim `defaultnsim'
    local seedm = cond(`seed' < 0, ".", "`seed'")
    if "`se'" == "" local se model
    local se = lower("`se'")
    if !inlist("`se'", "model", "robust") {
        di as err "se() must be model or robust"
        exit 198
    }
    if "`separation'" == "" local separation fail
    local separation = lower("`separation'")
    if !inlist("`separation'", "fail", "keep") {
        di as err "separation() must be fail or keep"
        exit 198
    }
    if "`name'" == "" local name _last
    sreturn local claim `claim'
    sreturn local sesoi `sesoi'
    sreturn local alpha `alpha'
    sreturn local alpharule `"`alpharule'"'
    sreturn local conf `conf'
    sreturn local nsim `nsim'
    sreturn local seed `seedm'
    sreturn local se `se'
    sreturn local separation `separation'
    sreturn local name `name'
    sreturn local disp = ("`display'" == "")
    sreturn local draws = ("`draws'" == "")
    sreturn local engine `engine'
    sreturn local effects `effects'
end

* post the result locals (set by the Mata layer in the caller's scope)
* usage: pape_post_res  (expects r_* locals via c_local trick: we instead
* call it as a subroutine that receives everything by option)

* ---------------------------------------------------------------------------
* power
* ---------------------------------------------------------------------------
program powerape_power, rclass
    syntax name(name=dname id="DGP name"), N(real) [*]
    pape_simopts, `options' defaultclaim(minimum) defaultnsim(1000)
    local n = floor(`n')
    pape_setengine, engine(`s(engine)') effects(`s(effects)')
    local cmd pape_cmd_power("`dname'", `n', "`s(claim)'", `s(sesoi)', `s(alpha)', `s(conf)', `"`s(alpharule)'"', `s(nsim)', `s(seed)', "`s(se)'", "`s(separation)'", `s(draws)', "`s(name)'", `s(disp)')
    pape_runmata `"`cmd'"' restore
    pape_post
    return add
end

program powerape_curve, rclass
    syntax name(name=dname id="DGP name"), N(numlist) [*]
    pape_simopts, `options' defaultclaim(minimum) defaultnsim(1000)
    pape_setengine, engine(`s(engine)') effects(`s(effects)')
    local cmd pape_cmd_curve("`dname'", "`n'", "`s(claim)'", `s(sesoi)', `s(alpha)', `s(conf)', `"`s(alpharule)'"', `s(nsim)', `s(seed)', "`s(se)'", "`s(separation)'", "`s(name)'", `s(disp)')
    pape_runmata `"`cmd'"' restore
    pape_post
    return add
end

program powerape_nreq, rclass
    syntax name(name=dname id="DGP name") [, POWer(real 0.90) NRANGE(numlist min=2 max=2) ///
        MAXiter(integer 5) noCONFirm NSIMConfirm(real -1) *]
    pape_simopts, `options' defaultclaim(minimum) defaultnsim(1500)
    if "`nrange'" == "" local nrange 30 2e6
    local n_lo : word 1 of `nrange'
    local n_hi : word 2 of `nrange'
    local conf_flag = ("`confirm'" == "")
    local nsim `s(nsim)'
    if `nsimconfirm' < 0 local nsimconfirm = 4 * `nsim'
    pape_setengine, engine(`s(engine)') effects(`s(effects)')
    local cmd pape_cmd_n("`dname'", `power', "`s(claim)'", `s(sesoi)', `s(alpha)', `s(conf)', `"`s(alpharule)'"', `nsim', `s(seed)', `n_lo', `n_hi', `maxiter', `conf_flag', `nsimconfirm', "`s(se)'", "`s(separation)'", "`s(name)'", `s(disp)')
    pape_runmata `"`cmd'"' restore
    pape_post
    return add
end

program powerape_mde, rclass
    syntax name(name=dname id="DGP name"), N(real) [POWer(real 0.80) MAINFocal(string) ///
        MAINModerator(string) MAXiter(integer 6) noCONFirm NSIMConfirm(real -1) DIRection(string) *]
    pape_simopts, `options' defaultclaim(detect) defaultnsim(1000)
    local n = floor(`n')
    if "`mainfocal'" == "" local mainfocal .
    if "`mainmoderator'" == "" local mainmoderator .
    local direction = lower("`direction'")
    if !inlist("`direction'", "", "positive", "negative") {
        di as err "direction() must be positive or negative"
        exit 198
    }
    local conf_flag = ("`confirm'" == "")
    local nsim `s(nsim)'
    if `nsimconfirm' < 0 local nsimconfirm = 4 * `nsim'
    pape_setengine, engine(`s(engine)') effects(`s(effects)')
    local cmd pape_cmd_mde("`dname'", `n', `power', "`s(claim)'", `s(sesoi)', `s(alpha)', `s(conf)', `"`s(alpharule)'"', `nsim', `s(seed)', `mainfocal', `mainmoderator', `maxiter', `conf_flag', `nsimconfirm', "`s(se)'", "`direction'", "`s(separation)'", "`s(name)'", `s(disp)')
    pape_runmata `"`cmd'"' restore
    pape_post
    return add
end

program powerape_robust, rclass
    syntax name(name=dname id="DGP name"), N(real) VARY(string) [PIN(string) MODE(string) ///
        POWer(real 0.80) GRIDpoints(integer 3) noNMAX NMAXPower(real 0.90) DIRection(string) *]
    pape_simopts, `options' defaultclaim(minimum) defaultnsim(800)
    local n = floor(`n')
    if "`pin'" == "" local pin ape
    local pin = lower("`pin'")
    if !inlist("`pin'", "ape", "coefficients") {
        di as err "pin() must be ape or coefficients"
        exit 198
    }
    if "`mode'" == "" local mode power
    local mode = lower("`mode'")
    if !inlist("`mode'", "power", "mde") {
        di as err "mode() must be power or mde"
        exit 198
    }
    local direction = lower("`direction'")
    * parse vary(name(values) name(values) ...)
    local rest = strtrim(`"`vary'"')
    local vnames ""
    local vvals ""
    while `"`rest'"' != "" {
        gettoken nm rest : rest, parse("(")
        local nm = lower(strtrim("`nm'"))
        gettoken vals rest : rest, match(par)
        if `"`par'"' != "(" {
            di as err "vary(): expected name(values) pairs"
            exit 198
        }
        local vnames `vnames' `nm'
        local vvals `"`vvals'`vals'|"'
        local rest = strtrim(`"`rest'"')
    }
    local nmaxflag = ("`nmax'" == "")
    pape_setengine, engine(`s(engine)') effects(`s(effects)')
    local cmd pape_cmd_robust("`dname'", `n', "`s(claim)'", `s(sesoi)', `s(alpha)', `s(conf)', `"`s(alpharule)'"', `s(nsim)', `s(seed)', "`vnames'", `"`vvals'"', "`pin'", "`mode'", `power', `gridpoints', `nmaxflag', `nmaxpower', "`s(se)'", "`direction'", "`s(separation)'", "`s(name)'", `s(disp)')
    pape_runmata `"`cmd'"' restore
    return matrix scenarios = pape_r_scenarios
    if "${pape_r_has_cand}" == "1" return matrix nmax_candidates = pape_r_nmaxcand
    foreach v in worst nmax nmax_power nmax_mcse nmax_scenario {
        capture return scalar `v' = scalar(pape_r_`v')
    }
    foreach v in fit engine effects {
        return local `v' `"${pape_r_`v'}"'
    }
    pape_clear_r
    return local vary "`vnames'"
    return local pin "`pin'"
    return local mode "`mode'"
end

program powerape_alpha, rclass
    syntax [name(name=resname id="result name")] [, COST(real 1) PRIOR(real 1) ERRor(string) ///
        CAP(real 0.05) SIZE NSIMSize(real -1) SEED(integer -1) Name(string) noDISplay ENGine(string) EFFects(string)]
    if "`resname'" == "" local resname _last
    if "`error'" == "" local error minimize
    local error = lower("`error'")
    if !inlist("`error'", "minimize", "balance") {
        di as err "error() must be minimize or balance"
        exit 198
    }
    local sizeflag = ("`size'" != "")
    local nsimsize = cond(`nsimsize' < 0, ".", "`nsimsize'")
    local seedm = cond(`seed' < 0, ".", "`seed'")
    if "`name'" == "" local name _lastalpha
    local disp = ("`display'" == "")
    pape_setengine, engine(`engine') effects(`effects')
    local cmd pape_cmd_alpha("`resname'", `cost', `prior', "`error'", `cap', `sizeflag', `nsimsize', `seedm', "`name'", `disp')
    pape_runmata `"`cmd'"' restore
    foreach v in alpha power mcse beta wcer power_default wcer_default range_lo range_hi ///
        capped cap alpha_uncapped wcer_uncapped size size_mcse size_truth {
        capture return scalar `v' = scalar(pape_r_`v')
    }
    return matrix grid = pape_r_grid
    pape_clear_r
    return local error "`error'"
    return local name "`name'"
end

program powerape_powerat, rclass
    syntax [name(name=resname id="result name")], ALPha(numlist)
    if "`resname'" == "" local resname _last
    mata: pape_cmd_powerat("`resname'", "`alpha'")
    return matrix table = pape_r_table
end

program powerape_frontier, rclass
    syntax [, COST(real 1) POWer(real 0.80) PRIOR(real 1) CLaim(string)]
    if "`claim'" == "" local claim detect
    local claim = lower("`claim'")
    if !inlist("`claim'", "detect", "minimum") {
        di as err "claim() must be detect or minimum"
        exit 198
    }
    mata: pape_cmd_frontier(`cost', `power', `prior', "`claim'")
    return scalar alpha = scalar(pape_r_alpha)
    pape_clear_r
end

program powerape_statement, rclass
    syntax [name(name=resname id="result name")]
    if "`resname'" == "" local resname _last
    mata: pape_cmd_statement("`resname'")
    return local statement `"${pape_r_statement}"'
    pape_clear_r
end

* ---------------------------------------------------------------------------
* post the r_* locals that the Mata layer sets in the CALLER's scope. Because
* locals are private to each program, the simulation commands call pape_post
* which reads them through c_local-free means: the Mata layer sets locals in
* the program that issued the mata: call (pape_runmata), so pape_runmata
* forwards them with c_local.
* ---------------------------------------------------------------------------
program pape_post, rclass
    foreach v in power mcse n nsim alpha conf sesoi target separated min_cell n_failed goal confirmed nsim_confirm mde is_margin {
        capture return scalar `v' = scalar(pape_r_`v')
    }
    foreach v in name kind claim estimand model se separation direction dgp alpharule engine effects fit {
        return local `v' `"${pape_r_`v'}"'
    }
    if "${pape_r_has_outcomes}" == "1" return matrix outcomes = pape_r_outcomes
    if "${pape_r_has_history}" == "1" return matrix history = pape_r_history
    if "${pape_r_has_results}" == "1" return matrix results = pape_r_results
    pape_clear_r
end
