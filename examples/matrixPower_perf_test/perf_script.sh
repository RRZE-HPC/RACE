############# Settings ###################

matrixFolder="/home/vault/unrz/unrz002h/matrix/MatrixPower"
folder="mtxPower_results"
powers="1 2 3 4 6 8 10"
threads="20"
nodes="2"
cacheSizes="15 20 25 35 45"
RCMs="0 1"

############ Don't change anything below ############################

./check-state.sh config.txt

mkdir -p $folder
mkdir -p $folder/raw_tradPower
mkdir -p $folder/raw_mtxPower

#get machine env
./machine-state.sh > $folder/machine-state.txt

#get matrix names
cd $matrixFolder
matrix_name=

while read -r i
do
    matrix_name=$matrix_name" "$i
done < <(find *)

echo $matrix_name

cd -


tradPower_file="${folder}/tradPower.csv"
printf "%20s, %10s, %5s, %5s, %5s\n" "# Matrix" "Perf" "Thread" "Power" "RCM" > ${tradPower_file}


function readResult
{
    cat $1 | grep "$2" | head -n 1 | cut -d "=" -f 2 | cut -d "G" -f1
}

for matrix in $matrix_name; do
    raw_file="${folder}/raw_tradPower/${matrix}.txt"
    for power in $powers; do
        for thread in $threads; do
            for RCM in $RCMs; do
                rcmFlag=""
                if [[ $RCM == "1" ]]; then
                    rcmFlag="-R"
                fi

                OMP_NUM_THREADS=$thread likwid-pin -c 0-$((thread-1)) ../tradPower \
                    -m ${matrixFolder}/${matrix} -c $thread -t 1 \
                    -i $power -v $rcmFlag > tmp.txt

                cat tmp.txt >> ${raw_file}
                perf=$(readResult "tmp.txt" "SpMV perf.")
                printf "%20s, %10.4f, %5d, %5d, %5d\n" ${matrix} ${perf} ${thread} ${power} ${RCM} >> ${tradPower_file}
            done
        done
    done
done


mtxPower_file="${folder}/mtxPower.csv"
printf "%20s, %10s, %5s, %5s, %5s, %5s, %5s\n" "# Matrix" "Perf" "Thread" "Power" "CS" "N" "RCM" > ${mtxPower_file}


for matrix in $matrix_name; do
    raw_file="${folder}/raw_mtxPower/${matrix}.txt"
    for power in $powers; do
        for thread in $threads; do
            for node in $nodes; do
                for cacheSize in $cacheSizes; do
                    cacheSizePerNode=$(echo $cacheSize/$node | bc -l)
                    for RCM in $RCMs; do
                        rcmFlag=""
                        if [[ $RCM == "1" ]]; then
                            rcmFlag="-R"
                        fi

                        OMP_NUM_THREADS=$thread likwid-pin -c 0-$((thread-1)) ../mtxPower \
                            -m ${matrixFolder}/${matrix} -c $thread -t 1 \
                            -i $power -n $node -s $cacheSizePerNode -v $rcmFlag > tmp.txt

                        cat tmp.txt >> ${raw_file}
                        perf=$(readResult "tmp.txt" "RACE power perf.")
                        printf "%20s, %10.4f, %5d, %5d, %5d, %5d, %5d\n" ${matrix} ${perf} ${thread} ${power} ${cacheSize} ${node} ${RCM} >> ${mtxPower_file}
                    done
                done
            done
        done
    done
done

rm -rf tmp.txt
