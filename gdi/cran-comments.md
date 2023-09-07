version 1.5.4:

## locally running R CMD check --as-cran
* using R version 4.3.1 (2023-06-16)
* using platform: x86_64-pc-linux-gnu (64-bit)
* R was compiled by
    gcc (Debian 10.2.1-6) 10.2.1 20210110
    GNU Fortran (Debian 10.2.1-6) 10.2.1 20210110
* running under: Debian GNU/Linux 11 (bullseye)
* using session charset: UTF-8
* using option ‘--as-cran’
* checking for file ‘gdi/DESCRIPTION’ ... OK
* this is package ‘gdi’ version ‘1.5.4’
* package encoding: UTF-8
* checking CRAN incoming feasibility ... [12s/16s] Note_to_CRAN_maintainers
Maintainer: ‘Darius Nau <dariusnau@gmx.at>’
* checking package namespace information ... OK
* checking package dependencies ... OK
* checking if this is a source package ... OK
* checking if there is a namespace ... OK
* checking for executable files ... OK
* checking for hidden files and directories ... OK
* checking for portable file names ... OK
* checking for sufficient/correct file permissions ... OK
* checking serialization versions ... OK
* checking whether package ‘gdi’ can be installed ... OK
* checking installed package size ... OK
* checking package directory ... OK
* checking for future file timestamps ... OK
* checking ‘build’ directory ... OK
* checking DESCRIPTION meta-information ... OK
* checking top-level files ... OK
* checking for left-over files ... OK
* checking index information ... OK
* checking package subdirectories ... OK
* checking R files for non-ASCII characters ... OK
* checking R files for syntax errors ... OK
* checking whether the package can be loaded ... OK
* checking whether the package can be loaded with stated dependencies ... OK
* checking whether the package can be unloaded cleanly ... OK
* checking whether the namespace can be loaded with stated dependencies ... OK
* checking whether the namespace can be unloaded cleanly ... OK
* checking loading without being on the library search path ... OK
* checking use of S3 registration ... OK
* checking dependencies in R code ... OK
* checking S3 generic/method consistency ... OK
* checking replacement functions ... OK
* checking foreign function calls ... OK
* checking R code for possible problems ... OK
* checking Rd files ... OK
* checking Rd metadata ... OK
* checking Rd line widths ... OK
* checking Rd cross-references ... OK
* checking for missing documentation entries ... OK
* checking for code/documentation mismatches ... OK
* checking Rd \usage sections ... OK
* checking Rd contents ... OK
* checking for unstated dependencies in examples ... OK
* checking installed files from ‘inst/doc’ ... OK
* checking files in ‘vignettes’ ... OK
* checking examples ... OK
* checking for unstated dependencies in ‘tests’ ... OK
* checking tests ... OK
  Running ‘testthat.R’
* checking for unstated dependencies in vignettes ... OK
* checking package vignettes in ‘inst/doc’ ... OK
* checking re-building of vignette outputs ... [16s/16s] OK
* checking PDF version of manual ... OK
* checking HTML version of manual ... OK
* checking for non-standard things in the check directory ... OK
* checking for detritus in the temp directory ... OK
* DONE
Status: OK

## devtools::check()

── R CMD check results ────────────────────────────────────────────────────────── gdi 1.5.4 ────
Duration: 54.3s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

## devtools::check_rhub()

─ gdi 1.5.4: NOTE

  Build ID:   gdi_1.5.4.tar.gz-50a62f5a404143e0a585178812e4728c
  Platform:   Windows Server 2022, R-devel, 64 bit
  Submitted:  3m 33.7s ago
  Build time: 3m 30.3s

❯ checking for non-standard things in the check directory ... NOTE
  Found the following files/directories:
    ''NULL''

❯ checking for detritus in the temp directory ... NOTE
  Found the following files/directories:
    'lastMiKTeXException'
       

0 errors ✔ | 0 warnings ✔ | 2 notes ✖


Build ID: 	gdi_1.5.4.tar.gz-fb8b7e0bcbee48c78e95a8f8dbc3d071
Platform: 	Fedora Linux, R-devel, clang, gfortran
Submitted: 	12 minutes 13.8 seconds ago
Build time: 	12 minutes 11.5 seconds
NOTES:

* checking HTML version of manual ... NOTE
Skipping checking HTML validation: no command 'tidy' found


Build ID: 	gdi_1.5.4.tar.gz-bdc365a72276473bb9ea66255d6f81d0
Platform: 	Ubuntu Linux 20.04.1 LTS, R-release, GCC
Submitted: 	12 minutes 47.7 seconds ago
Build time: 	12 minutes 45.4 seconds
NOTES:

* checking HTML version of manual ... NOTE
Skipping checking HTML validation: no command 'tidy' found
