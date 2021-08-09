############# Settings ###################

matrixFolder="/home/vault/unrz/unrz002h/matrix/MatrixPower_1"
folder="mtxDetails_1"
thread=64
############ Don't change anything below ############################

mkdir -p $folder

#get matrix names
cd $matrixFolder
matrix_name=

while read -r i
do
    matrix_name=$matrix_name" "$i
done < <(find *)

echo $matrix_name

cd -


tradPower_file="${folder}/details.csv"
printf "%20s, %18s, %18s\n" "# Matrix" "Nrows" "Nnz" > ${tradPower_file}

for matrix in $matrix_name; do
    OMP_NUM_THREADS=$thread likwid-pin -c 0-$((thread-1)) ../matStat \
        -m ${matrixFolder}/${matrix} -c $thread -t 1 > tmp.txt

    nrows_w_space=$(grep "Nrows =" tmp.txt | cut -d"=" -f2)
    nrows=$(echo ${nrows_w_space})
    nnz_w_space=$(grep "NNZ =" tmp.txt | cut -d"=" -f2)
    nnz=$(echo ${nnz_w_space})


    printf "%20s, %18d, %18d\n" ${matrix} ${nrows} ${nnz} >> ${tradPower_file}
done

rm -rf tmp.txt
