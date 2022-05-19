#!/bin/bash

cp template_power_all.tex power.tex
cd powerPlots
str=""
for matrix in *.tex; do
    str="${str}\\\input{$matrix} "
done
../substitute.sh "Input_files" "${str}" ../power.tex
mv ../power.tex .
cd -
