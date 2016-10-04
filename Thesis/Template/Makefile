TEXFILES=$(wildcard *.tex)
PDFFILES=$(patsubst %.tex,%.pdf,$(TEXFILES))
DVIFILES=$(patsubst %.tex,%.dvi,$(TEXFILES))
AUXFILES=$(patsubst %.tex,%.aux,$(TEXFILES))
LOGFILES=$(patsubst %.tex,%.log,$(TEXFILES))
OUTFILES=$(patsubst %.tex,%.out,$(TEXFILES))

all: $(PDFFILES) $(DVIFILES)

%.pdf: %.tex clean
	pdflatex $*
	pdflatex $*

%.dvi: %.tex clean
	latex $*
	latex $*
	
clean:
	rm -f $(AUXFILES) $(LOGFILES) $(OUTFILES)

distclean: clean
	rm -f $(PDFFILES) $(DVIFILES)

.PHONY: all cleanall
