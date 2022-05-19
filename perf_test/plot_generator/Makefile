pdf:
	pdflatex symmspmv.tex
	cd powerPlots && pdflatex power.tex cd -
	make clean_tmp

read:
	evince *.pdf &

clean_tmp:
	rm -f *.log *.ps *.aux *.out *.dvi *.bbl *.blg *.alg *.glo *.glg *.acr *.gls *.ist *.acn
	cd powerPlots && rm -f *.log *.ps *.aux *.out *.dvi *.bbl *.blg *.alg *.glo *.glg *.acr *.gls *.ist *.acn && cd -

clean:
	make clean_tmp
	rm -f *.pdf
	rm -f powerPlots/*.pdf
