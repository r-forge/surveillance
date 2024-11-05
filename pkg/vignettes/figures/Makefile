### create static figures with extra dependencies

all: HMM.pdf fPeriods2.pdf fPeriods3.pdf

## PStricks-based figure for surveillance.Rnw
HMM.pdf: HMM.tex
	latex -interaction=batchmode $< && dvips -q ${<:.tex=.dvi} && \
		ps2pdf -q -sOutputFile=- ${<:.tex=.ps} | pdfcrop - $@
	@rm -f $(addprefix $(basename $<).,aux dvi log ps)

## fig:fPlot in monitoringCounts.Rnw
fPeriods2.pdf fPeriods3.pdf &: fPeriods.R
	R --vanilla -s -f "$<"

.PHONY: all
