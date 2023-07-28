version 1.4.1:
>>fixed a gap in the previous versions concerning single-channel images

version 1.4.0:
>>The primary addition to this version is a plotting function, plot_sil(), for visualizing silhouettes used with the main functions of the package, and two functions for calculating the horizontal and vertical position of the center of mass (COM), COM() and hCOM(). Estimating COM positions is often of importance for drawing inferences about posture and locomotion in extinct organism, e.g. determining whether an animal was quadrupedal or bipedal.

## locally running R CMD check --as-cran
* checking for file ‘gdi/DESCRIPTION’ ... OK
* preparing ‘gdi’:
* checking DESCRIPTION meta-information ... OK
* installing the package to build vignettes
* creating vignettes ... OK
* checking for LF line-endings in source and make files and shell scripts
* checking for empty or unneeded directories
* building ‘gdi_1.4.1.tar.gz’

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
* this is package ‘gdi’ version ‘1.4.1’
* package encoding: UTF-8
* checking CRAN incoming feasibility ... [11s/17s] Note_to_CRAN_maintainers
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
* checking tests ...
  Running ‘testthat.R’
 OK
* checking for unstated dependencies in vignettes ... OK
* checking package vignettes in ‘inst/doc’ ... OK
* checking re-building of vignette outputs ... OK
* checking PDF version of manual ... OK
* checking HTML version of manual ... OK
* checking for non-standard things in the check directory ... OK
* checking for detritus in the temp directory ... OK
* DONE

Status: OK


## devtools::check()
─  using R version 4.3.1 (2023-06-16)
─  using platform: x86_64-pc-linux-gnu (64-bit)
─  R was compiled by
       gcc (Debian 10.2.1-6) 10.2.1 20210110
       GNU Fortran (Debian 10.2.1-6) 10.2.1 20210110
─  running under: Debian GNU/Linux 11 (bullseye)
─  using session charset: UTF-8
─  using options ‘--no-manual --as-cran’
✔  checking for file ‘gdi/DESCRIPTION’ ...
─  this is package ‘gdi’ version ‘1.4.1’
─  package encoding: UTF-8
✔  checking package namespace information ...
✔  checking package dependencies (1.8s)
✔  checking if this is a source package ...
✔  checking if there is a namespace
✔  checking for executable files ...
✔  checking for hidden files and directories
✔  checking for portable file names
✔  checking for sufficient/correct file permissions
✔  checking serialization versions
✔  checking whether package ‘gdi’ can be installed (3.5s)
✔  checking installed package size
✔  checking package directory
✔  checking for future file timestamps
✔  checking ‘build’ directory
✔  checking DESCRIPTION meta-information ...
✔  checking top-level files
✔  checking for left-over files
✔  checking index information (406ms)
✔  checking package subdirectories (477ms)
✔  checking R files for non-ASCII characters ...
✔  checking R files for syntax errors ...
✔  checking whether the package can be loaded (439ms)
✔  checking whether the package can be loaded with stated dependencies ...
✔  checking whether the package can be unloaded cleanly ...
✔  checking whether the namespace can be loaded with stated dependencies ...
✔  checking whether the namespace can be unloaded cleanly (396ms)
✔  checking loading without being on the library search path (453ms)
✔  checking dependencies in R code (740ms)
✔  checking S3 generic/method consistency (550ms)
✔  checking replacement functions ...
✔  checking foreign function calls (391ms)
✔  checking R code for possible problems (6.9s)
✔  checking Rd files (573ms)
✔  checking Rd metadata ...
✔  checking Rd line widths ...
✔  checking Rd cross-references ...
✔  checking for missing documentation entries ...
✔  checking for code/documentation mismatches (963ms)
✔  checking Rd \usage sections (928ms)
✔  checking Rd contents ...
✔  checking for unstated dependencies in examples ...
✔  checking installed files from ‘inst/doc’ ...
✔  checking files in ‘vignettes’ ...
✔  checking examples (5.6s)
✔  checking for unstated dependencies in ‘tests’ ...
─  checking tests (398ms)
✔  Running ‘testthat.R’ (5.2s)
✔  checking for unstated dependencies in vignettes (444ms)
✔  checking package vignettes in ‘inst/doc’ ...
─  checking re-building of vignette outputs ... [10s/10s] OK (10.4s)
✔  checking for non-standard things in the check directory
✔  checking for detritus in the temp directory
   
   
── R CMD check results ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────── gdi 1.4.1 ────
Duration: 45.6s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔


## devtools::check_rhub()
── gdi 1.4.1: NOTE

  Build ID:   gdi_1.4.1.tar.gz-2f91f9004a924da1899a85e28d8e01a1
  Platform:   Windows Server 2022, R-devel, 64 bit
  Submitted:  3m 20.8s ago
  Build time: 3m 14.2s

❯ checking for non-standard things in the check directory ... NOTE
  Found the following files/directories:
    ''NULL''

❯ checking for detritus in the temp directory ... NOTE
  Found the following files/directories:
    'lastMiKTeXException'

0 errors ✔ | 0 warnings ✔ | 2 notes ✖

── gdi 1.4.1: IN-PROGRESS

  Build ID:   gdi_1.4.1.tar.gz-6774370508fa43e38b390e4b63291e79
  Platform:   Ubuntu Linux 20.04.1 LTS, R-release, GCC
  Submitted:  3m 20.9s ago


── gdi 1.4.1: IN-PROGRESS

  Build ID:   gdi_1.4.1.tar.gz-dc2b820e59c641789cd8cbac313c78e1
  Platform:   Fedora Linux, R-devel, clang, gfortran
  Submitted:  3m 20.9s ago

  
## devtools::check_win_devel()

