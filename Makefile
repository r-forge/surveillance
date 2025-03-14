################################################################################
## Build and check the package
## CAVE: Watch the list in .Rbuildignore
##
## Authors: Sebastian Meyer, Michael Hoehle
## $Date$
##
## History:
##    8 Jan 2008 (MH): initial version
##   29 Mar 2012 (SM): Added making of R/sysdata.rda from sysdata/sysdata.R
##    4 Apr 2012 (SM): R execution command as a variable
##   11 Jun 2012 (SM): Restored and cleaned up Makefile
##   18 Jun 2012 (SM): fixed make manual bug (must first remove old manual.pdf)
##   19 Jun 2012 (SM): run check on built package instead of source directory
##   27 Jun 2012 (SM): added checkUsage recipe (R package codetools)
##    9 Aug 2012 (SM): added --timings for R CMD check
##   27 Aug 2012 (SM): install from built .tar.gz (which includes vignettes)
##   14 Sep 2012 (SM): automatic date in DESCRIPTION file
##   15 Nov 2012 (SM): removed automatic date in DESCRIPTION file since R-Forge
##                     now includes date and revision in its built packages
##   10 Dec 2012 (SM): add message about warnings in examples in "check"
##    6 Feb 2014 (SM): set _R_CHECK_FORCE_SUGGESTS_ to FALSE (for INLA package)
##   25 Aug 2014 (SM): use compact-vignettes=both, i.e., including ghostscript
##    7 Jun 2015 (SM): added rule to create the NEWS.html as on CRAN
##   12 Jun 2015 (SM): added rule to check with allExamples
##   17 Mar 2016 (SM): check-allExamples with --run-dontrun and --extra-arch
##   21 Jun 2017 (SM): account for R code with roxygen documentation
##   12 Jul 2017 (SM): "quick" vs. CRAN-versions of build and check rules
##   13 Sep 2018 (SM): drop roxygen (no longer supports latin1 packages)
##   02 Nov 2018 (SM): add pkgdown recipe
##               ... : (see svn log)
################################################################################

## Define variable for R which enables the use of alternatives,
## e.g., 'make check R=R-devel'
R := R
export LANGUAGE := en

## package version
VERSION := $(shell $R --vanilla -s -e 'cat(read.dcf("pkg/DESCRIPTION", fields="Version"))')

## build the package
build:
	$R CMD build --no-resave-data --compact-vignettes=both pkg
build-noVignettes:
	$R CMD build --no-resave-data --no-build-vignettes pkg

## Save internal datasets from pkg/sysdata/ into pkg/R/sysdata.rda
sysdata:
	cd pkg/sysdata; $R CMD BATCH --vanilla --no-timing sysdata.R
	mv pkg/sysdata/sysdata.rda pkg/R/sysdata.rda


## auxiliary functions ("canned recipes") for check rules
define CHECK_REPORT_TIMINGS_SCRIPT
timings <- read.table(file.path('surveillance.Rcheck','surveillance-Ex.timings'), header=TRUE, row.names='name')
timings <- timings[order(timings[['elapsed']], decreasing=TRUE),]
cat("Runtimes of examples with elapsed time > 2 seconds:\n\n")
subset(timings, elapsed > 2)
endef
export CHECK_REPORT_TIMINGS_SCRIPT
define check-report-timings
echo "$${CHECK_REPORT_TIMINGS_SCRIPT}" | $R -s --vanilla
endef

define check-report-warnings-in-examples
(cd surveillance.Rcheck; \
nwarn=`grep -cP '^Warning(?!.*(gpc.poly|pit.default|k = 40|newly enabled))' surveillance-Ex.Rout`; \
if [ $$nwarn -gt 0 ]; then echo "\n\tWARNING: $$nwarn" \
	"warning(s) thrown when running examples,\n" \
	"\t         see file surveillance.Rcheck/surveillance-Ex.Rout\n"; fi)
endef

## "quick" check
check: build-noVignettes
	_R_CHECK_FORCE_SUGGESTS_=FALSE _R_CHECK_COMPACT_DATA_=FALSE \
	_R_CHECK_PKG_SIZES_=FALSE _R_CHECK_DOC_SIZES_=FALSE \
	$R CMD check --no-manual --ignore-vignettes --check-subdirs=no surveillance_${VERSION}.tar.gz

## standard --as-cran check
check-cran: export _R_CHECK_PKG_SIZES_THRESHOLD_ := 8.5
check-cran: build
	$R CMD check --as-cran --timings surveillance_${VERSION}.tar.gz
## further option: --use-gct (for better detection of memory bugs/segfaults)
	@$(check-report-timings)
	@$(check-report-warnings-in-examples)

## check with "allExamples" (donttest), --run-demo and --run-dontrun
## also use --extra-arch to only do runtime tests (no R and Rd code checking)
check-allExamples: build-noVignettes
	_R_SURVEILLANCE_ALL_EXAMPLES_=TRUE $R CMD check --ignore-vignettes --timings --run-demo --run-dontrun --extra-arch --output=/tmp surveillance_${VERSION}.tar.gz
	@cd /tmp; $(check-report-timings); $(check-report-warnings-in-examples)


install: build
	$R CMD INSTALL surveillance_${VERSION}.tar.gz

checkUsage: install
	echo "library('surveillance'); library('codetools'); \
	checkUsagePackage('surveillance', suppressFundefMismatch=FALSE, \
	    suppressLocalUnused=TRUE, suppressNoLocalFun=TRUE, skipWith=TRUE, \
	    suppressUndefined=FALSE, suppressPartialMatchArgs=FALSE)" \
	| $R -s --no-save --no-restore

spelling:
	cd pkg; \
	  codespell -q 35 -S '*~' -L 'ans,parm,hist' NEWS.md README.md R man vignettes

## we need to run Rd2pdf inside pkg such that \packageTitle finds DESCRIPTION
manuals:
	$R CMD Rd2pdf --batch --force --output=manual.pdf pkg
	cd pkg; \
	  echo "tools::pkg2HTML(dir = '.', out = '../manual.html', toc_entry = 'name')" \
	  | $R -s --no-save --no-restore
	xdg-open manual.html

NEWS.html: pkg/NEWS.md
README.html: pkg/README.md
%.html: pkg/%.md
	pandoc -f markdown -t html -s --css=https://CRAN.R-project.org/web/CRAN_web.css -o "$@" "$<"

www:
	cd pkg; $R --no-echo --no-save --no-restore -e \
	  "pkgdown::build_site(examples = FALSE, lazy = TRUE, devel = TRUE); summary(warnings())"

www-clean:
	cd pkg; $R --no-echo --no-save --no-restore -e \
	  "pkgdown::clean_site(); pkgdown::build_site(examples = FALSE, new_process = FALSE); summary(warnings())"

www-applications: www/applications_EE.csv
	cd www; $R --no-echo --no-save --no-restore -e \
	  "rmarkdown::render('applications_EE.Rmd')"

clean:
	make -C pkg/demo clean
	cd pkg/src; rm -f *.o *.so *.dll symbols.rds
	make -C pkg/vignettes clean
	rm -f pkg/*/.Rhistory

.PHONY: build build-noVignettes sysdata check check-cran check-allExamples install checkUsage spelling manuals www clean
