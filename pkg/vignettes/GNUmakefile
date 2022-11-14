## This is a convenient Makefile for quick vignette compilation in Emacs
## and for maintainer-level checks of the vignette cache. It requires GNU make.
## Because of this and because it does nothing special to build
## the vignettes, it is listed in the top-level .Rbuildignore file
## and called GNUmakefile so that it is not used by R CMD build.
## In addition it defines a rule to create the PStricks-based HMM figure
## included in surveillance.Rnw
##
## Copyright (C) 2016,2022  Sebastian Meyer <seb.meyer@fau.de>
## This Makefile can be redistributed and/or modified under the terms
## of the GNU General Public License, version 2 or later.

R := R
RNWFILES := $(wildcard *.Rnw)
TEXFILESFROMRNW := ${RNWFILES:.Rnw=.tex}
VIGNETTES := ${RNWFILES:.Rnw=.pdf}
CACHEFILES := $(wildcard $(patsubst %.Rnw,%-cache.RData,${RNWFILES}))

.PHONY: all
all: ${VIGNETTES}

%.pdf: %.tex references.bib
	${R} --vanilla -s -e "tools::texi2pdf('$<', clean = TRUE)"

monitoringCounts.pdf: monitoringCounts.bib  # additional dependency

.SECONDARY: ${TEXFILESFROMRNW}  # make shall not remove intermediate tex files
%.tex: %.Rnw
	LANGUAGE=en ${R} --vanilla CMD Sweave "$<"

.PHONY: clean cleanall
clean:
	rm -rf *.aux *.bbl *.blg *.log *.out *-R.bib plots/ Rplots.pdf $(addprefix *-cache.RData,.new .diff) *.Rout
cleanall: clean
	rm -rf ${VIGNETTES} ${TEXFILESFROMRNW} ${RNWFILES:.Rnw=.R}


### create surveillance-hmm.pdf

#surveillance.pdf: surveillance-hmm.pdf  # additional dependency
surveillance-hmm.pdf: surveillance-hmm.tex
	latex -interaction=batchmode $< && dvips -q ${<:.tex=.dvi} && \
		ps2pdf -q -sOutputFile=- ${<:.tex=.ps} | pdfcrop - $@
	rm -f $(addprefix $(basename $<).,aux dvi log ps)


### check running tangled vignette code because knitr tampers with R CMD check,
### see <https://github.com/yihui/knitr/pull/2036>

.PHONY: check-tangle
check-tangle:
	@rm -f ${RNWFILES:.Rnw=.Rout} ${RNWFILES:.Rnw=.R}  # force re-purl()
	$(MAKE) ${RNWFILES:.Rnw=.Rout}
%.Rout: %.R
	LANGUAGE=en ${R} CMD BATCH --vanilla "$<" && rm -f Rplots.pdf || { tail -13 "$@"; exit 1; }
%.R: %.Rnw
	LANGUAGE=en ${R} --vanilla -s -e 'COMPUTE <- FALSE; tools::buildVignette("$<", weave=FALSE, tangle=TRUE)'


### rebuild cache and check if it differs from the old one

.PHONY: check-cache
check-cache: $(addsuffix .diff,${CACHEFILES})
	ls -l $^

$(addsuffix .diff,${CACHEFILES}): %.diff: % %.new
	@echo; echo CHECKING $<.new AGAINST $< ...
	@${R} --vanilla -s -e "\
old <- attach('$<'); \
new <- attach('$<.new', warn.conflicts = FALSE); \
all.equal.ahg <- function(target,current,...) {current$$time <- target$$time; NextMethod()}; \
equal <- all.equal(old, new); \
if (isTRUE(equal)) message('OK: $< needs no update') else { \
message('CAVE: $< MAY NEED TO BE UPDATED\n-> $<.diff'); print(equal) }" > "$@"

$(addsuffix .new,${CACHEFILES}): %.new: %
	@echo; echo REBUILDING $< AS $@ ...
	mv "$<" "$<.bak"
	-${MAKE} -B ${<:-cache.RData=.tex}
	-mv "$<" "$@"
	mv "$<.bak" "$<"
