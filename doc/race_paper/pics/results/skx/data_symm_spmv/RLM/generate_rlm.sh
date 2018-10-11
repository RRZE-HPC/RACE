#!/bin/bash


NNZ_TABLE=table_nnz.txt
NNZ_COL="1 2"
NNZ_DELIM=" "

SPMV_TABLE=../MKL/result.txt
SPMV_COL="2 3"
SPMV_DELIM="|"

ALPHA_TABLE=SpMV_RCM/result.txt
ALPHA_COL="2 4"
ALPHA_DELIM="|"

OUT_FILE=../RLM/result.txt

bw_copy=104
bw_load=115

#prepare nnz table
col1=$(echo $NNZ_COL | cut -d" " -f1) 
a=$(cut -d"$NNZ_DELIM" -f $col1 $NNZ_TABLE)
echo $a | tr " " "\n" > tmp1.txt
col2=$(echo $NNZ_COL | cut -d" " -f2) 
b=$(cut -d"$NNZ_DELIM" -f $col2 $NNZ_TABLE)
echo $b | tr " " "\n" > tmp2.txt
paste -d" " tmp1.txt tmp2.txt > tmp3.txt
sort -k1 -t"," tmp3.txt > sorted_nnz.txt

#prepare spmv table, with removed header
col1=$(echo $SPMV_COL | cut -d" " -f1) 
a=$(cut -d"$SPMV_DELIM" -f $col1 $SPMV_TABLE | tail -n +2)
echo $a | tr " " "\n" > tmp1.txt
col2=$(echo $SPMV_COL | cut -d" " -f2) 
b=$(cut -d"$SPMV_DELIM" -f $col2 $SPMV_TABLE | tail -n +2)
echo $b | tr " " "\n" > tmp2.txt
paste -d" " tmp1.txt tmp2.txt > tmp3.txt
sort -k1 -t"," tmp3.txt > sorted_spmv.txt

#prepare measured alpha table, with removed header
col1=$(echo $ALPHA_COL | cut -d" " -f1) 
a=$(cut -d"$ALPHA_DELIM" -f $col1 $ALPHA_TABLE | tail -n +2)
echo $a | tr " " "\n" > tmp1.txt
col2=$(echo $ALPHA_COL | cut -d" " -f2)
b=$(cut -d"$ALPHA_DELIM" -f $col2 $ALPHA_TABLE | tail -n +2)
echo $b | tr " " "\n" > tmp2.txt
paste -d" " tmp1.txt tmp2.txt > tmp3.txt
sort -k1 -t"," tmp3.txt > sorted_measured_alpha.txt

rm tmp1.txt tmp2.txt tmp3.txt

echo "Id | Matrix | NNZR | SymmSpMV measured alpha copy | SymmSpMV opt alpha copy | SymmSpMV SpMV alpha copy | SymmSpMV measured alpha ld | SymmSpMV opt alpha ld | SymmSpMV SpMV alpha ld" > $OUT_FILE
printf "Id,\t Matrix,\t\t\t NNZR,\t Measured alpha\n" > alphas.txt
julia rlm_model.jl $bw_copy $bw_load sorted_nnz.txt sorted_spmv.txt sorted_measured_alpha.txt tmp.txt
cat tmp.txt >> $OUT_FILE
cat alpha.txt >> alphas.txt
rm sorted_nnz.txt sorted_spmv.txt tmp.txt alpha.txt

