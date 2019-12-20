filename=race

#pdf: ps
#	ps2pdf ${filename}.ps

#pdf-print: ps
#	ps2pdf -dColorConversionStrategy=/LeaveColorUnchanged -dPDFSETTINGS=/printer ${filename}.ps

#text: html
#	html2text -width 100 -style pretty ${filename}/${filename}.html | sed -n '/./,$$p' | head -n-2 >${filename}.txt

#html:
#	@#latex2html -split +0 -info "" -no_navigation ${filename}
#	htlatex ${filename}

#ps:	dvi
#	dvips -t letter ${filename}.dvi

pdf:
	pdflatex ${filename}
	makeglossaries ${filename}
	bibtex ${filename}||true
	pdflatex ${filename}
	pdflatex ${filename}

read:
	evince ${filename}.pdf &

#aread:
#	acroread ${filename}.pdf

clean:
	rm -f ${filename}.log ${filename}.pdf ${filename}.ps ${filename}.aux ${filename}.out ${filename}.dvi ${filename}.bbl ${filename}.blg
