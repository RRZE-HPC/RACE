#!/bin/bash

#sort according to row
./permute_w_rows.sh

curr_folder=$PWD

for arch in ivy skx; do
	for kernel in symm_kacz symm_spmv; do 
		echo "ARCH = $arch"
		cd $arch/data_"$kernel"/plot_generator/

		folder_names=
		while read -r i
		do
			folder_names=$folder_names" "$i
		done < <(find * -type d)


		for folder in $folder_names; do
			#find script and config file
			cd $folder
			#script_file=$(find * -name "*.sh")
			config_file=$(find * -name "*.txt")
			newConfig_file=$config_file".perm"
			#use permuted result file
			cat $config_file | sed -e "s@result.txt@result_permuted.txt@g">$newConfig_file
			script_file=$curr_folder/generate_plot.sh
			echo "executing" $script_file $newConfig_file
			$script_file $newConfig_file
			#remove new config file and permuted results
			rm $newConfig_file
			cd -
		done

		cd $curr_folder
	done
done

cd skx/corner_cases_scaling
./generate_plots.sh

cd -

./rm_permute_files.sh

