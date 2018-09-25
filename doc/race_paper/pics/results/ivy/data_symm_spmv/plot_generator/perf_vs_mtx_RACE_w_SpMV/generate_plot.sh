#!/bin/bash

col_mtx_name=`cat $1 | grep COL_MTX_NAMES  | awk '{print $3}'`
template=`cat $1 | grep TEMPLATE  | awk '{print $3}'`
generated=`cat $1 | grep GENERATED  | awk '{print $3}'`
col_perf_1=`cat $1 | grep COL_PERF_1 | awk '{ for(i=3;i<=NF;i++) {print $i} }'`
col_perf_2=`cat $1 | grep COL_PERF_2 | awk '{ for(i=3;i<=NF;i++) {print $i} }'`
file_op=`cat $1 | grep OP | awk '{print $3}'`
perf_files=`cat $1 | grep PERF_FILES | awk '{ for(i=3;i<=NF;i++) {print $i} }'`
colors=`cat $1 | grep COLORS | awk '{ for(i=3;i<=NF;i++) {print $i} }'`
title=`cat $1 | grep TITLE | awk '{ for(i=3;i<=NF;i++) {print $i} }'`
legends=`cat $1 | grep LEGEND | awk '{ for(i=3;i<=NF;i++) {print $i} }'`
ylabel=`cat $1 | grep YLABEL | awk '{ for(i=3;i<=NF;i++) {print $i} }'`
type=`cat $1 | grep TYPE | awk '{print $3}'`
yscale=`cat $1 | grep Y_SCALE | awk '{print $3}'`
yscale_file=`cat $1 | grep Y_SCALE_FILE | awk '{print $3}'`
yscale_col=`cat $1 | grep Y_SCALE_COL | awk '{print $3}'`
yscale_op=`cat $1 | grep Y_SCALE_OP | awk '{print $3}'`

mark="default triangle square diamond otimes"

mtx_name_file=$(echo $perf_files | cut -d " " -f  1)
cat $mtx_name_file | tail -n +2 >temp.txt

sort -k$col_mtx_name temp.txt>sorted_file.txt
rm temp.txt
mtx_name_file="sorted_file.txt"


cut -f$col_mtx_name $mtx_name_file > col0.txt
sed -i -- "s/|//g" col0.txt
sed -i -- "s/_/-/g" col0.txt


counter=1

if [ "$yscale_file" != "" ]; then
	cut -f$yscale_col $yscale_file | grep   -Eo '[0-9]*\.?[0-9]+'>"scale_file.txt"
fi

for file in $perf_files ; do
	#sort perf files
	sort -t"|" -k$col_mtx_name $file>temp.txt
	file="temp.txt"
	curr_col=$(echo $col_perf_1 | cut -d " " -f  $counter)
	cut -d "|" -f$curr_col $file | grep   -Eo '[0-9]*\.?[0-9]+'>"col_1_$counter.txt"

	if [ "$col_perf_2" != "" ]; then
		curr_col_2=$(echo $col_perf_2 | cut -d " " -f  $counter)
		cut -d "|" -f$curr_col_2 $file | grep   -Eo '[0-9]*\.?[0-9]+'>"col_2_$counter.txt"
		paste -d$file_op col_1_$counter.txt col_2_$counter.txt | bc -l > "col$counter.txt"
		rm "col_1_$counter.txt" "col_2_$counter.txt"
	else
		mv col_1_$counter.txt col$counter.txt
	fi

	if [ "$yscale_file" != "" ]; then
		cp "col$counter.txt" "col$counter_tmp.txt"
		paste -d* col$counter_tmp.txt scale_file.txt | bc -l >col$counter.txt
		rm col$counter_tmp.txt 
	fi

	rm temp.txt
	let counter=counter+1
done

tot_matrices=$(wc -l <col0.txt) 
#generate latex file
#correct x-ticks
cp $template $generated
#generate xticks
xticks=1
for id in `seq 2 $tot_matrices`; do
    xticks=$xticks", $id"
done

counter=1
mtx_name=
while IFS='' read -r line || [[ -n "$line" ]]; do
 	base_name=`basename $line .mtx`
	if [ $counter == 1 ] 
	then
		mtx_name=$base_name
	else
	    	mtx_name=$mtx_name", $base_name"
	fi

let counter=counter+1
done < "col0.txt"

echo $mtx_name

sed -i -- "s/xtick=.*/xtick={$xticks},/g" $generated
sed -i -- "s/xticklabels=.*/xticklabels={$mtx_name},/g" $generated

touch temp_coo.txt
col_counter=1
#generate coordinates
for file in $perf_files ; do
	coo=
	counter=1
	while IFS='' read -r line || [[ -n "$line" ]]; do
		if [ $counter == 1 ] 
		then
			coo="($counter,$line)"
		else
		    	coo=$coo" ($counter,$line)"
		fi

	let counter=counter+1
	done < "col$col_counter.txt"

	curr_color=$(echo $colors | cut -d " " -f  $col_counter)

	#sed -i -- "s/plot coordinates.*/plot coordinates{$coo};,/g" $generated

	y_scale_code=""

	if [ "$yscale" != "" ]; then
		if [ "$yscale_file" = "" ]; then
			y_scale_code=", y filter/.code={\pgfmathparse{\pgfmathresult*$yscale}\pgfmathresult}"
		fi
	fi

	if [ $type == "SCATTER" ]; then
		let xshift=($col_counter-1)*30
		curr_mark=$(echo $mark | cut -d " " -f  $col_counter)
		if [ "$curr_mark" == "default" ] 
		then
			curr_mark=
		fi
		echo $curr_mark
		echo "\\\addplot[xshift=$xshift, mark="$curr_mark*", only marks, mark size=10pt, fill=$curr_color,draw=$curr_color $y_scale_code] plot coordinates{$coo};">>temp_coo.txt
	elif [ $type == "LINE" ]; then
		let xshift=($col_counter-1)*30
		curr_mark=$(echo $mark | cut -d " " -f  $col_counter)
		if [ "$curr_mark" == "default" ]
		then
			curr_mark=
		fi
		echo $curr_mark
		echo "\\\addplot[mark="$curr_mark*", mark size=10pt, mark options={$curr_color}, draw=$curr_color $y_scale_code] plot coordinates{$coo};">>temp_coo.txt

	else
		echo "\\\addplot[ybar,bar width=0.7cm,fill=$curr_color,draw=$curr_color $y_scale_code] plot coordinates{$coo};">>temp_coo.txt
	fi

	let col_counter=col_counter+1
done

	coo=$(cat temp_coo.txt)

	#find line number for substitution
	line_num=$(grep -rne "addplot cmd" $template | cut -f 1 | grep   -Eo '[0-9]*\.?[0-9]+')
	echo  $line_num
	awk -v n=$line_num -v s="$coo" 'NR == n {print s} {print}' $generated > temp.tex
	
	#old_plot_cmd="addplot cmd"
	#cat $generated | eval $(echo "sed 's@$old_plot_cmd@$coo@g'") >temp.tex
	rm $generated
	mv temp.tex $generated

	#sed -i -- "s/addplot cmd/a $coo/g" $generated
	
	#add legend
	counter=1
	legend=
	for currLegend in $legends; do
		if [ $counter == 1 ] 
		then
			legend=$currLegend
		else
		    	legend=$legend", $currLegend"
		fi
		let counter=counter+1
	done
	echo $legend
	
	sed -i -- "s/legend comes here/$legend/g" $generated

	echo $title
	#add title
	old_title="title comes here"
	cat $generated | eval $(echo "sed 's/$old_title/$title/g'") > temp_title.txt
	rm $generated
	mv temp_title.txt $generated

	echo $ylabel
	old_ylabel="ylabel comes here"
	cat $generated | eval $(echo "sed 's#$old_ylabel#$ylabel#g'") > temp_title.txt
	rm $generated
	mv temp_title.txt $generated

	bar_commands="%add_bar_commands"
	if [ $type != "LINE" ]; then
		bar_commands="ybar=2*\\\pgflinewidth,bar width=14pt,"
	fi
	old_bar_commands="%add_bar_commands"
	cat $generated | eval $(echo "sed 's@$old_bar_commands@$bar_commands@g'") > temp_title.txt
	rm $generated
	mv temp_title.txt $generated
	#sed -i -- "s/title comes here/$title/g" $generated

	legend_commands="%spl_legend_code"
	if [ $type != "LINE" ]; then
		legend_commands="legend image code/.code={\\\draw[#1, draw=none] (0cm,-0.2cm) rectangle (0.6cm,0.2cm);}," 
	fi
	old_legend_commands="%spl_legend_code"
	cat $generated | eval $(echo "sed 's@$old_legend_commands@$legend_commands@g'") > temp_title.txt
	rm $generated
	mv temp_title.txt $generated

counter=1
for file in $perf_files ; do
	rm "col$counter.txt"
	let counter=counter+1
done

rm sorted_file.txt
rm temp_coo.txt
rm col0.txt

#make
rm *.aux
rm *.bbl
rm *.blg
rm *.log
