* build.do -- compile the powerape Mata sources into lpowerape.mlib
* Repository use: from the folder that holds build.do and ado/ (a clone of
*   https://github.com/jespernwulff/powerape-stata), run  do build.do
*   It compiles the sources in ado/ and writes the library there.
* Installed use (after -net install powerape- and -net get powerape-, which copies
*   this file to the current directory): do build.do
*   -net install- puts the .mata sources on the adopath (PLUS); build.do finds them
*   with -findfile- and replaces lpowerape.mlib where the installed copy lives (or
*   writes it to the current directory if none is installed). Then run
*   mata: mata mlib index. Needed on Stata 16-18: the shipped library was compiled
*   by StataNow 19.5 and a library compiled by a newer Stata may not load in an
*   older one.
* Stata 16 or later; lmbuild is an official command ([M-3] lmbuild).
version 16
clear all
local here "`c(pwd)'"
local root "`here'"
local sources powerape_core powerape_structs powerape_engine powerape_panel powerape_iv ///
    powerape_stata_engine powerape_search powerape_cmd
capture confirm file "`root'/ado/powerape_core.mata"
if !_rc {
    local src "`root'/ado"
    local out "`root'/ado"
}
else {
    * not a copy of the repository: the sources are on the adopath (net install)
    capture findfile powerape_core.mata
    if _rc {
        di as err "build.do: powerape_core.mata not found in ado/ below the current folder or on the adopath"
        exit 601
    }
    local fn `"`r(fn)'"'
    local src = substr(`"`fn'"', 1, strlen(`"`fn'"') - strlen("powerape_core.mata") - 1)
    capture findfile lpowerape.mlib
    if _rc local out "`here'"
    else {
        local fn `"`r(fn)'"'
        local out = substr(`"`fn'"', 1, strlen(`"`fn'"') - strlen("lpowerape.mlib") - 1)
    }
}
di as txt "build.do: sources in `src'"
di as txt "build.do: library goes to `out'"
mata: mata clear
foreach f of local sources {
    do "`src'/`f'.mata"
}
capture erase "`out'/lpowerape.mlib"
lmbuild lpowerape, dir("`out'") replace
mata: mata mlib index
di as txt "lpowerape.mlib built in `out'"
cd "`here'"
