export LC_NUMERIC="en_US.UTF-8"
############# Settings ###################

matrixFolder="/home/vault/unrz/unrz002h/matrix/MatrixPower"
folder="mtxPower_results_w_likwid"
threads="20"
paramFileRACE="paramRACE.csv"
paramFileTRAD="paramTRAD.csv"
trad="on"

############ Don't change anything below ############################

./check-state.sh config.txt

mkdir -p $folder
mkdir -p $folder/raw_tradPower
mkdir -p $folder/raw_mtxPower

#get machine env
./machine-state.sh > $folder/machine-state.txt

matrix_name=$(cut -d ";" -f1 ${paramFileRACE} | tail -n +2)
function readResult
{
    cat $1 | grep "$2" | head -n 1 | cut -d "=" -f 2 | cut -d "G" -f1
}

function readDataVol
{
    string="$1"
    file=$2
    selector=$3

    iter=$(cat ${file} | grep "Num iterations =" | cut -d"=" -f2)
    data_vol=$(grep "${string}" ${file} | head -n ${selector} | tail -n 1 | cut -d"," -f 2)

    data_vol_per_iter=$(echo "${data_vol}*1000/${iter}" | bc -l) #in MB
    echo ${data_vol_per_iter}
}

powers=$(cut -d ";" -f2 ${paramFileTRAD} | tail -n +2)
RCMs=$(cut -d ";" -f3 ${paramFileTRAD} | tail -n +2)

if [[ ${trad} == "on" ]]; then
    tradPower_file="${folder}/tradPower.csv"
    printf "%20s, %12s, %12s, %10s, %8s, %8s, %5s, %5s, %5s, %8s, %8s, %8s, %8s, %8s, %8s\n" "#Matrix" "NROWS" "NNZ" "Perf" "Clock" "Uncore" "Thread" "Power" "RCM" "MEM_rd" "MEM_wr" "L3_rd" "L3_wr" "L2_rd" "L2_wr" > ${tradPower_file}
    ctr=1
    for matrix in $matrix_name; do
        raw_file="${folder}/raw_tradPower/${matrix}.txt"
        for thread in $threads; do
            power=$(echo ${powers} | cut -d" " -f${ctr})
            RCM=$(echo ${RCMs} | cut -d" " -f${ctr})
            rcmFlag=""
            if [[ $RCM == "1" ]]; then
                rcmFlag="-R"
            fi

            OMP_NUM_THREADS=$thread likwid-pin -c 0-$((thread-1)) ../tradPower \
                -m ${matrixFolder}/${matrix} -c $thread -t 1 \
                -i $power -v $rcmFlag &> tmp.txt

            cat tmp.txt >> ${raw_file}
            perf=$(readResult "tmp.txt" "SpMV perf.")
            nrows=$(cat tmp.txt | grep "Nrows =" | cut -d"=" -f2 | cut -d"," -f1)
            nnz=$(cat tmp.txt | grep "Nrows =" | cut -d"=" -f3 | cut -d"," -f1)

            selector=1 #which region
            OMP_NUM_THREADS=$thread likwid-perfctr -m -g MEM -O --stats -C 0-$((thread-1)) ../tradPower \
                -m ${matrixFolder}/${matrix} -c $thread -t 1 \
                -i $power -v $rcmFlag > tmp.txt
            cat tmp.txt >> ${raw_file}
            mem_read_trad=$(readDataVol "Memory read data volume \[GBytes\] STAT" tmp.txt ${selector})
            mem_write_trad=$(readDataVol "Memory write data volume \[GBytes\] STAT" tmp.txt ${selector})

            OMP_NUM_THREADS=$thread likwid-perfctr -m -g L3 -O --stats -C 0-$((thread-1)) ../tradPower \
                -m ${matrixFolder}/${matrix} -c $thread -t 1 \
                -i $power -v $rcmFlag > tmp.txt
            cat tmp.txt >> ${raw_file}
            l3_read_trad=$(readDataVol "L3 load data volume \[GBytes\] STAT" tmp.txt ${selector})
            l3_write_trad=$(readDataVol "L3 evict data volume \[GBytes\] STAT" tmp.txt ${selector})

            OMP_NUM_THREADS=$thread likwid-perfctr -m -g L2 -O --stats -C 0-$((thread-1)) ../tradPower \
                -m ${matrixFolder}/${matrix} -c $thread -t 1 \
                -i $power -v $rcmFlag > tmp.txt
            cat tmp.txt >> ${raw_file}
            l2_read_trad=$(readDataVol "L2D load data volume \[GBytes\] STAT" tmp.txt ${selector})
            l2_write_trad=$(readDataVol "L2D evict data volume \[GBytes\] STAT" tmp.txt ${selector})

            OMP_NUM_THREADS=$thread likwid-perfctr -m -g CLOCK -O --stats -C 0-$((thread-1)) ../tradPower \
                -m ${matrixFolder}/${matrix} -c $thread -t 1 \
                -i $power -v $rcmFlag > tmp.txt
            cat tmp.txt >> ${raw_file}
            clock=$(grep "^Clock \[MHz\] STAT" tmp.txt | head -n ${selector} | tail -n 1 | cut -d"," -f 5)
            uncore=$(grep "^Uncore Clock \[MHz\] STAT" tmp.txt | head -n ${selector} | tail -n 1 | cut -d"," -f 2)

            printf "%20s, %12d, %12d, %10.4f, %8.3f, %8.3f, %5d, %5d, %5d, %8.3f, %8.3f, %8.3f, %8.3f, %8.3f, %8.3f\n" ${matrix} ${nrows} ${nnz} ${perf} ${clock} ${uncore} ${thread} ${power} ${RCM} ${mem_read_trad} ${mem_write_trad} ${l3_read_trad} ${l3_write_trad} ${l2_read_trad} ${l2_write_trad}>> ${tradPower_file}
        done
        let ctr=${ctr}+1
    done
fi


mtxPower_file="${folder}/mtxPower.csv"
printf "%20s, %12s, %12s, %10s, %10s, %8s, %8s, %8s, %8s, %5s, %5s, %5s, %5s, %50s, %5s, %8s, %8s, %8s, %8s, %8s, %8s, %8s, %8s, %8s, %8s, %8s, %8s\n" "# Matrix" "NROWS" "NNZ" "Perf" "TRAD_perf" "ClockRACE" "UncoreRACE" "ClockTrad" "UncoreTrad" "Thread" "Power" "CS" "N" "SF" "RCM" "MEM_rd_race" "MEM_wr_race" "L3_rd_race" "L3_wr_race" "L2_rd_race" "L2_wr_race" "MEM_rd_trad" "MEM_wr_trad" "L3_rd_trad" "L3_wr_trad" "L2_rd_trad" "L2_wr_trad" > ${mtxPower_file}

powers=$(cut -d ";" -f2 ${paramFileRACE} | tail -n +2)
cacheSizes=$(cut -d ";" -f3 ${paramFileRACE} | tail -n +2)
RCMs=$(cut -d ";" -f4 ${paramFileRACE} | tail -n +2)
stages=$(cut -d ";" -f5 ${paramFileRACE} | tail -n +2)

ctr=1
for matrix in $matrix_name; do
    raw_file="${folder}/raw_mtxPower/${matrix}.txt"
    for thread in $threads; do
        power=$(echo ${powers} | cut -d" " -f${ctr})
        RCM=$(echo ${RCMs} | cut -d" " -f${ctr})
        cacheSize=$(echo ${cacheSizes} | cut -d" " -f${ctr})
        stage=$(echo ${stages} | cut -d" " -f${ctr})
        node=1
        safetyFactor=""

        for((idx=0; idx<${stage}; ++idx)); do
            safetyFactor="${safetyFactor}1,"
        done
        safetyFactor="${safetyFactor}100000"

        cacheSizePerNode=$(echo $cacheSize/$node | bc -l)
        rcmFlag=""
        if [[ $RCM == "1" ]]; then
            rcmFlag="-R"
        fi

        RACE_CACHE_VIOLATION_CUTOFF=${safetyFactor} OMP_NUM_THREADS=$thread likwid-pin -c 0-$((thread-1)) ../mtxPower \
            -m ${matrixFolder}/${matrix} -c $thread -t 1 \
            -i $power -n $node -s $cacheSizePerNode -v $rcmFlag &> tmp.txt

        cat tmp.txt >> ${raw_file}
        spmv_perf=$(readResult "tmp.txt" "SpMV power perf.")
        perf=$(readResult "tmp.txt" "RACE power perf.")
        nrows=$(cat tmp.txt | grep "Nrows =" | cut -d"=" -f2 | cut -d"," -f1)
        nnz=$(cat tmp.txt | grep "Nrows =" | cut -d"=" -f3 | cut -d"," -f1)

        RACE_CACHE_VIOLATION_CUTOFF=${safetyFactor} OMP_NUM_THREADS=$thread likwid-perfctr -m -g MEM -O --stats -C 0-$((thread-1)) ../mtxPower \
            -m ${matrixFolder}/${matrix} -c $thread -t 1 \
            -i $power -n $node -s $cacheSizePerNode -v $rcmFlag > tmp.txt
        cat tmp.txt >> ${raw_file}
        selector=1 #which region
        mem_read_trad=$(readDataVol "Memory read data volume \[GBytes\] STAT" tmp.txt ${selector})
        mem_write_trad=$(readDataVol "Memory write data volume \[GBytes\] STAT" tmp.txt ${selector})
        selector=2 #which region
        mem_read_race=$(readDataVol "Memory read data volume \[GBytes\] STAT" tmp.txt ${selector})
        mem_write_race=$(readDataVol "Memory write data volume \[GBytes\] STAT" tmp.txt ${selector})


        RACE_CACHE_VIOLATION_CUTOFF=${safetyFactor} OMP_NUM_THREADS=$thread likwid-perfctr -m -g L3 -O --stats -C 0-$((thread-1)) ../mtxPower \
            -m ${matrixFolder}/${matrix} -c $thread -t 1 \
            -i $power -n $node -s $cacheSizePerNode -v $rcmFlag > tmp.txt
        cat tmp.txt >> ${raw_file}
        selector=1
        l3_read_trad=$(readDataVol "L3 load data volume \[GBytes\] STAT" tmp.txt ${selector})
        l3_write_trad=$(readDataVol "L3 evict data volume \[GBytes\] STAT" tmp.txt ${selector})
        selector=2
        l3_read_race=$(readDataVol "L3 load data volume \[GBytes\] STAT" tmp.txt ${selector})
        l3_write_race=$(readDataVol "L3 evict data volume \[GBytes\] STAT" tmp.txt ${selector})


        RACE_CACHE_VIOLATION_CUTOFF=${safetyFactor} OMP_NUM_THREADS=$thread likwid-perfctr -m -g L2 -O --stats -C 0-$((thread-1)) ../mtxPower \
            -m ${matrixFolder}/${matrix} -c $thread -t 1 \
            -i $power -n $node -s $cacheSizePerNode -v $rcmFlag > tmp.txt
        cat tmp.txt >> ${raw_file}
        selector=1
        l2_read_trad=$(readDataVol "L2D load data volume \[GBytes\] STAT" tmp.txt ${selector})
        l2_write_trad=$(readDataVol "L2D evict data volume \[GBytes\] STAT" tmp.txt ${selector})
        selector=2
        l2_read_race=$(readDataVol "L2D load data volume \[GBytes\] STAT" tmp.txt ${selector})
        l2_write_race=$(readDataVol "L2D evict data volume \[GBytes\] STAT" tmp.txt ${selector})


        RACE_CACHE_VIOLATION_CUTOFF=${safetyFactor} OMP_NUM_THREADS=$thread likwid-perfctr -m -g CLOCK -O --stats -C 0-$((thread-1)) ../mtxPower \
            -m ${matrixFolder}/${matrix} -c $thread -t 1 \
            -i $power -n $node -s $cacheSizePerNode -v $rcmFlag > tmp.txt
        cat tmp.txt >> ${raw_file}
        selector=1
        clock_trad=$(grep "^Clock \[MHz\] STAT" tmp.txt | head -n ${selector} | tail -n 1 | cut -d"," -f 5)
        uncore_trad=$(grep "^Uncore Clock \[MHz\] STAT" tmp.txt | head -n ${selector} | tail -n 1 | cut -d"," -f 2)
        selector=2
        clock_race=$(grep "^Clock \[MHz\] STAT" tmp.txt | head -n ${selector} | tail -n 1 | cut -d"," -f 5)
        uncore_race=$(grep "^Uncore Clock \[MHz\] STAT" tmp.txt | head -n ${selector} | tail -n 1 | cut -d"," -f 2)

        safetyFactor_wo_comma=$(echo ${safetyFactor} | sed -e "s@,@;@g")
        printf "%20s, %12d, %12d, %10.4f, %10.4f, %8.3f, %8.3f, %8.3f, %8.3f, %5d, %5d, %5.2f, %5d, %50s, %5d, %8.3f, %8.3f, %8.3f, %8.3f, %8.3f, %8.3f, %8.3f, %8.3f, %8.3f, %8.3f, %8.3f, %8.3f\n" ${matrix} ${nrows} ${nnz} ${perf} ${spmv_perf} ${clock_race} ${uncore_race} ${clock_trad} ${uncore_trad} ${thread} ${power} ${cacheSizePerNode} ${node} ${safetyFactor_wo_comma} ${RCM} ${mem_read_race} ${mem_write_race} ${l3_read_race} ${l3_write_race} ${l2_read_race} ${l2_write_race} ${mem_read_trad} ${mem_write_trad} ${l3_read_trad} ${l3_write_trad} ${l2_read_trad} ${l2_write_trad} >> ${mtxPower_file}
    done
    let ctr=${ctr}+1
done

rm -rf tmp.txt
