FIGS := sa4 sa3-population
FIGS += sa3-deltas sa3-deltas-map
FIGS += influenza-case-distribution influenza-classify-median
FIGS += sa3-min-err-map risk-ratio-cmp

ANALYSIS := data/model-analysis/model-errors-sa3.ssv

PDFS := $(FIGS:%=%.pdf)

all: $(ANALYSIS) $(PDFS)

clean:
	@rm -f $(PDFS) data/model-analysis/*.ssv

.PHONY: all clean

$(ANALYSIS): analyse-output.R data/model-output
	@./analyse-output.R

sa4.pdf: data/mixing/sa3-populations.ssv map.R
	@./map.R --discrete --height 5 --name SA4 --title SA4 --annotate $< $@

sa3-population.pdf: data/mixing/sa3-populations.ssv map.R
	@./map.R --continuous --width 6 --height 5 --name pcnt_of_popn --percentage $< $@

sa3-deltas.pdf: data/model-analysis/model-errors-sa3.ssv plot-sa3-deltas.R
	@./plot-sa3-deltas.R $< $@

sa3-deltas-map.pdf: data/model-analysis/sa3-deltas.ssv map.R
	@./map.R --discrete --name delta --facet label --width 10.75 --height 9 $< $@

influenza-case-distribution.pdf: data/influenza/influenza-case-distribution.ssv plot-flu-distr.R
	@./plot-flu-distr.R $< $@

influenza-classify-median.pdf:data/influenza/influenza-median-delta.ssv map.R
	@./map.R --discrete --name category --height 4.25 $< $@

risk-ratio-cmp.pdf: data/model-analysis/model-errors-sa3.ssv data/influenza/influenza-notif-ratio.ssv plot-risk-ratio-cmp.R
	@./plot-risk-ratio-cmp.R

sa3-min-err-map.pdf: data/model-analysis/sa3-min-err.ssv map.R
	@./map.R --diverging --name err --pad 3 --facet label --width 9.75 --height 9 $< $@
