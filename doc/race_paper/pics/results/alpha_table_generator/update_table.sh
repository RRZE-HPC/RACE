#!/bin/bash

julia update_table.jl
template_file=table_template.tex

#remove " signs
sed -i 's@"@@g' table.txt

#now permute table
./permute_w_rows.sh

#escape _ sign
sed -i 's@_@\\_@g' table_permuted.txt

#count number of total matrices
lines=$(cut -f 1 -d "&" table_permuted.txt | wc -l)

#make reordered id
for i in `seq 1 $lines`; do
	echo "{"$i"}">>key.txt
done

cut -f 2- -d "&" table_permuted.txt > rest.txt

paste -d"&" key.txt rest.txt > table_permuted.txt

rm key.txt rest.txt

#substitute in template
data=$(cat table_permuted.txt)
line_num=$(grep -rne "#TABLE_DATA#" $template_file | cut -f 1 | grep   -Eo '[0-9]*\.?[0-9]+')
awk -v n=$line_num -v s="$data" 'NR == n {print s} {print}' $template_file > table.tex


