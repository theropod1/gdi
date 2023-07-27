version 1.4.0:
>>The primary addition to this version is a plotting function, plot.sil(), for visualizing silhouettes used with the main functions of the package, and two functions for calculating the horizontal and vertical position of the center of mass (COM), COM() and hCOM(). Estimating COM positions is often of importance for drawing inferences about posture and locomotion in extinct organism, e.g. determining whether an animal was quadrupedal or bipedal.

## locally running R CMD check --as-cran
* checking for file ‘gdi/DESCRIPTION’ ... OK
* preparing ‘gdi’:
* checking DESCRIPTION meta-information ... OK
* installing the package to build vignettes
* creating vignettes ... OK
* checking for LF line-endings in source and make files and shell scripts
* checking for empty or unneeded directories
* building ‘gdi_1.4.0.tar.gz’

* using log directory ‘/home/suirad/Documents/research_projects/gdi_package/gdi.Rcheck’
* using R version 4.3.1 (2023-06-16)
* using platform: x86_64-pc-linux-gnu (64-bit)
* R was compiled by
    gcc (Debian 10.2.1-6) 10.2.1 20210110
    GNU Fortran (Debian 10.2.1-6) 10.2.1 20210110
* running under: Debian GNU/Linux 11 (bullseye)
* using session charset: UTF-8
* using option ‘--as-cran’
* checking for file ‘gdi/DESCRIPTION’ ... OK
* this is package ‘gdi’ version ‘1.4.0’
* package encoding: UTF-8
* checking CRAN incoming feasibility ... [10s/15s] Note_to_CRAN_maintainers
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
* checking examples ... [12s/12s] OK
* checking for unstated dependencies in ‘tests’ ... OK
* checking tests ...
  Running ‘testthat.R’
 OK
* checking for unstated dependencies in vignettes ... OK
* checking package vignettes in ‘inst/doc’ ... OK
* checking re-building of vignette outputs ... [18s/18s] OK
* checking PDF version of manual ... OK
* checking HTML version of manual ... OK
* checking for non-standard things in the check directory ... OK
* checking for detritus in the temp directory ... OK
* DONE

Status: OK


## devtools::check_rhub()
── gdi 1.4.0: NOTE

  Build ID:   gdi_1.4.0.tar.gz-7a1ff44b2ff24d53b1ad303675bdccb9
  Platform:   Windows Server 2022, R-devel, 64 bit
  Submitted:  3m 35.3s ago
  Build time: 3m 32.1s

❯ checking for non-standard things in the check directory ... NOTE
  

❯ checking for detritus in the temp directory ... NOTE
  Found the following files/directories:
    'lastMiKTeXException'

0 errors ✔ | 0 warnings ✔ | 2 notes ✖

── gdi 1.4.0: IN-PROGRESS

  Build ID:   gdi_1.4.0.tar.gz-da2c5541077b49e9889dc4522e6f34f7
  Platform:   Ubuntu Linux 20.04.1 LTS, R-release, GCC
  Submitted:  3m 35.4s ago


── gdi 1.4.0: IN-PROGRESS

  Build ID:   gdi_1.4.0.tar.gz-1d71a669d6d84598bec916b692b9eb9c
  Platform:   Fedora Linux, R-devel, clang, gfortran
  Submitted:  3m 35.4s ago

## devtools::check_win_devel()

* using log directory 'd:/RCompile/CRANguest/R-devel/gdi.Rcheck'
* using R Under development (unstable) (2023-07-26 r84754 ucrt)
* using platform: x86_64-w64-mingw32
* R was compiled by
    gcc.exe (GCC) 12.2.0
    GNU Fortran (GCC) 12.2.0
* running under: Windows Server 2022 x64 (build 20348)
* using session charset: UTF-8
* checking for file 'gdi/DESCRIPTION' ... OK
* this is package 'gdi' version '1.4.0'
* package encoding: UTF-8
* checking CRAN incoming feasibility ... Note_to_CRAN_maintainers
Maintainer: 'Darius Nau <dariusnau@gmx.at>'
* checking package namespace information ... OK
* checking package dependencies ... OK
* checking if this is a source package ... OK
* checking if there is a namespace ... OK
* checking for hidden files and directories ... OK
* checking for portable file names ... OK
* checking serialization versions ... OK
* checking whether package 'gdi' can be installed ... OK
* checking installed package size ... OK
* checking package directory ... OK
* checking for future file timestamps ... OK
* checking 'build' directory ... OK
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
* checking S3 generic/method consistency ... NOTE
Mismatches for apparent methods not registered:
plot:
  function(x, ...)
plot.sil:
  function(sil, flip, add, ...)
See section 'Registering S3 methods' in the 'Writing R Extensions'
manual.
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
* checking installed files from 'inst/doc' ... OK
* checking files in 'vignettes' ... OK
* checking examples ... OK
* checking for unstated dependencies in 'tests' ... OK
* checking tests ... OK
  Running 'testthat.R'
* checking for unstated dependencies in vignettes ... OK
* checking package vignettes in 'inst/doc' ... OK
* checking re-building of vignette outputs ... OK
* checking PDF version of manual ... [13s] OK
* checking HTML version of manual ... OK
* DONE
Status: 1 NOTE
