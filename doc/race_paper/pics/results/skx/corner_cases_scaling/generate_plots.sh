#!/bin/bash

MATRICES="crankseg_1 Graphene-4096 inline_1 parabolic_fem"
RESULTS="scaling_results"
OUT="plots"

OUT_NORCM="$OUT/NORCM"
OUT_RCM="$OUT/RCM"

mkdir -p plots_tmp
mkdir -p $OUT_NORCM
mkdir -p $OUT_RCM

for mat in $MATRICES; do
	for ext in NORCM RCM; do
		cp plot_template.tex plots_tmp/$mat.tex
		./substitute.sh "FOLDER" "$PWD/$RESULTS/$ext"  plots_tmp/$mat.tex
		./substitute.sh "MATRIX_NAME" "$mat" plots_tmp/$mat.tex
		cd plots_tmp
		echo "compiling $mat"
		pdflatex $mat.tex >> compile.txt
		echo "finished compiling $mat"
		cd -
		cp plots_tmp/$mat.pdf $OUT/$ext/.
	done
done

rm -r plots_tmp

