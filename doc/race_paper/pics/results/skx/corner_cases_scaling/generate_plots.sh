#!/bin/bash

MATRICES="crankseg_1 Graphene-4096 inline_1 parabolic_fem"
RLM_copy="-1 25.33 32.61  21.88"
RLM_load="-1 28.01 36.06 24.19"
RESULTS="scaling_results"
OUT="plots"

OUT_NORCM="$OUT/NORCM"
OUT_RCM="$OUT/RCM"

mkdir -p plots_tmp
mkdir -p $OUT_NORCM
mkdir -p $OUT_RCM

ctr=1

for mat in $MATRICES; do
	for ext in NORCM RCM; do

		rlm_copy=$(echo $RLM_copy | cut -f $ctr -d " ")
		rlm_load=$(echo $RLM_load | cut -f $ctr -d " ")
		template="plot_template.tex" 
		echo $rlm_copy
		if [ "$rlm_copy" != "-1" ] ; then
			template="plot_template_with_rlm.tex"
		fi
		cp $template plots_tmp/$mat.tex
		./substitute.sh "FOLDER" "$PWD/$RESULTS/$ext"  plots_tmp/$mat.tex
		./substitute.sh "MATRIX_NAME" "$mat" plots_tmp/$mat.tex
		if [ "$rlm_copy" != "-1" ] ; then
			./substitute.sh "RLM_COPY" "$rlm_copy" plots_tmp/$mat.tex
			./substitute.sh "RLM_LOAD" "$rlm_load" plots_tmp/$mat.tex
		fi
		cd plots_tmp
		echo "compiling $mat"
		pdflatex $mat.tex >> compile.txt
		echo "finished compiling $mat"
		cd -
		cp plots_tmp/$mat.pdf $OUT/$ext/.
	done
	let ctr=$ctr+1
done

rm -r plots_tmp

