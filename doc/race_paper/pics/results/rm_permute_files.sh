#!/bin/bash

curr_folder=$PWD

for arch in ivy skx; do
	for kernel in symm_kacz symm_spmv; do 
		for method in ABMC MC MKL RACE RSB; do
			echo "ARCH = $arch"
			permResultFile=$arch/data_"$kernel"/$method/result_permuted.txt
			rm -f $permResultFile
		done
	done
done


