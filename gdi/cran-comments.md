version 1.2.3:
>>minor fix

version 1.2.2:
>>Compared to version 1.1.1, v 1.2.2 supplies some minor fixes and updates and three additional convenience functions that help to automatically determine the correct background colour settings and allow the generation of coordinates for plotting and visually assessing superelliptical body cross-sections of various exponents.

## locally running R CMD check --as-cran
* using R version 4.3.0 (2023-04-21)
* using platform: x86_64-pc-linux-gnu (64-bit)
* R was compiled by
    gcc (Debian 10.2.1-6) 10.2.1 20210110
    GNU Fortran (Debian 10.2.1-6) 10.2.1 20210110
* running under: Debian GNU/Linux 11 (bullseye)
* using session charset: UTF-8
* using option ‘--as-cran’
* checking for file ‘gdi/DESCRIPTION’ ... OK
* this is package ‘gdi’ version ‘1.2.2’
* package encoding: UTF-8
* checking CRAN incoming feasibility ... [10s/19s] Note_to_CRAN_maintainers
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
* checking re-building of vignette outputs ... OK
* checking PDF version of manual ... OK
* checking HTML version of manual ... NOTE
Skipping checking HTML validation: no command 'tidy' found
* checking for non-standard things in the check directory ... OK
* checking for detritus in the temp directory ... OK
* DONE
Status: 1 NOTE


## devtools::check_rhub()
── gdi 1.2.2: NOTE

  Build ID:   gdi_1.2.2.tar.gz-b3984da3f4fd4abb882c8f1befe3a771
  Platform:   Windows Server 2022, R-devel, 64 bit
  Submitted:  3m 3.9s ago
  Build time: 2m 48.2s

❯ checking for non-standard things in the check directory ... NOTE
  Found the following files/directories:
    ''NULL''

❯ checking for detritus in the temp directory ... NOTE
  Found the following files/directories:
    'lastMiKTeXException'

0 errors ✔ | 0 warnings ✔ | 2 notes ✖


## devtools::check_win_devel()

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Darius Nau <dariusnau@gmx.at>'

Checking DOIs failed with message:
object 'ind' not found
[…]
Status: 1 NOTE
