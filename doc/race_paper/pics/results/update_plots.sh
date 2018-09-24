#!/bin/bash

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
			script_file=$(find * -name "*.sh")
			config_file=$(find * -name "*.txt")
			echo "executing" $script_file $config_file
			./$script_file $config_file
			cd -
		done

		cd $curr_folder
	done
done

