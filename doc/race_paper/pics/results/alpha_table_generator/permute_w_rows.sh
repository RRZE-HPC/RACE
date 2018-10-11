#!/bin/bash

curr_folder=$PWD

sort_col=3
mtxDetailFile=matrix_details.txt

resultFile=table.txt
permResultFile=table_permuted.txt
cp $resultFile $permResultFile
#echo "KEY" > tmp.tmp
echo $(cut -d"|" -f$sort_col $mtxDetailFile | grep [0-9])|tr " " "\n">tmp.tmp
paste -d"|" tmp.tmp $resultFile > tmp_file.tmp
head -n 1 tmp_file.tmp > tmp_file_2.tmp
(tail -n +2 tmp_file.tmp | sort -s -k1 -n ) >> tmp_file_2.tmp
#now remove first element
cut -d"|" -f2- tmp_file_2.tmp > $permResultFile
	
rm tmp_file.tmp tmp_file_2.tmp tmp.tmp


