export LC_NUMERIC="en_US.UTF-8"

#### This script does not distinguish between load and store traffic

############# Settings ###################
configFile=$1
matrixFolder=$(cat ${configFile} | grep "matrixFolder" | cut -d"=" -f2)
folder=$(cat ${configFile} | grep "folder" | cut -d"=" -f2)
threads=$(cat ${configFile} | grep "threads" | cut -d"=" -f2)
paramFileRACE=$(cat ${configFile} | grep "paramFileRACE" | cut -d"=" -f2)
paramFileTRAD=$(cat ${configFile} | grep "paramFileTRAD" | cut -d"=" -f2)
trad="on"

############ Don't change anything below ############################

./check-state.sh config.txt

mkdir -p $folder
mkdir -p $folder/raw_tradPower
mkdir -p $folder/raw_mtxPower

#get machine env
./machine-state.sh > $folder/machine-state_old.txt
MachineState/machinestate.py -o $folder/machine-state.json

matrix_name=$(cut -d "," -f1 ${paramFileRACE} | tail -n +2)
function readResult
{
    cat $1 | grep "$2" | head -n 1 | cut -d "=" -f 2 | cut -d "G" -f1
}

function readMetric
{
    string="$1"
    file=$2
    selector=$3

    data_vol=$(grep "${string}" ${file} | head -n ${selector} | tail -n 1 | cut -d"," -f 2)
    echo ${data_vol}
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

powers=$(cut -d "," -f2 ${paramFileTRAD} | tail -n +2)
RCMs=$(cut -d "," -f3 ${paramFileTRAD} | tail -n +2)

if [[ ${trad} == "on" ]]; then
    tradPower_file="${folder}/tradPower.csv"
    printf "%20s, %12s, %12s, %10s, %8s, %5s, %5s, %5s, %8s, %8s, %8s, %8s, %8s, %8s, %10s\n" "#Matrix" "NROWS" "NNZ" "Perf" "Clock" "Thread" "Power" "RCM" "MEM" "MEM_bw" "L3" "L3_bw" "L2" "L2_bw" "Total_Pre-time" > ${tradPower_file}
    ctr=1
    for matrix in $matrix_name; do
        raw_file="${folder}/raw_tradPower/${matrix}.txt"
        for thread in $threads; do
            power=$(echo ${powers} | cut -d" " -f${ctr})
            RCM_w_space=$(echo ${RCMs} | cut -d" " -f${ctr})
            RCM=$(echo $RCM_w_space)
            echo "$RCM"
            rcmFlag=""
            if [[ $RCM == "1" ]]; then
                rcmFlag="-R"
            fi

            echo "RCM = $rcmFlag"

            OMP_NUM_THREADS=$thread OMP_PROC_BIND=close OMP_PLACES=cores numactl -m 0 ../tradPower \
                -m ${matrixFolder}/${matrix} -c $thread -t 1 \
                -i $power -v $rcmFlag &> tmp.txt

            cat tmp.txt >> ${raw_file}
            perf=$(readResult "tmp.txt" "SpMV perf.")
            nrows=$(cat tmp.txt | grep "Nrows =" | cut -d"=" -f2 | cut -d"," -f1)
            nnz=$(cat tmp.txt | grep "Nrows =" | cut -d"=" -f3 | cut -d"," -f1)
            preTime=$(cat tmp.txt | grep "Total pre-processing time =" | cut -d"=" -f2 | cut -d"s" -f1)

            selector=1 #which region
            OMP_NUM_THREADS=$thread numactl -m 0 likwid-perfctr -m -g MEM -O -V 3 --stats -C 0-$((thread-1)) ../tradPower \
                -m ${matrixFolder}/${matrix} -c $thread -t 1 \
                -i $power -v $rcmFlag > tmp.txt
            cat tmp.txt >> ${raw_file}
            mem_data_trad=$(readDataVol "Memory data volume \[GBytes\] STAT" tmp.txt ${selector})
            mem_bw_trad=$(readMetric "Memory bandwidth \[MBytes/s\] STAT" tmp.txt ${selector})

            OMP_NUM_THREADS=$thread numactl -m 0 likwid-perfctr -m -g L3 -O -V 3 --stats -C 0-$((thread-1)) ../tradPower \
                -m ${matrixFolder}/${matrix} -c $thread -t 1 \
                -i $power -v $rcmFlag > tmp.txt
            cat tmp.txt >> ${raw_file}
            l3_data_trad=$(readDataVol "L3 data volume \[GBytes\] STAT" tmp.txt ${selector})
            l3_bw_trad=$(readMetric "L3 bandwidth \[MBytes/s\] STAT" tmp.txt ${selector})


            OMP_NUM_THREADS=$thread numactl -m 0 likwid-perfctr -m -g L2 -O -V 3 --stats -C 0-$((thread-1)) ../tradPower \
                -m ${matrixFolder}/${matrix} -c $thread -t 1 \
                -i $power -v $rcmFlag > tmp.txt
            cat tmp.txt >> ${raw_file}
            l2_data_trad=$(readDataVol "L2 data volume \[GBytes\] STAT" tmp.txt ${selector})
            l2_bw_trad=$(readMetric "L2 bandwidth \[MBytes/s\] STAT" tmp.txt ${selector})


            OMP_NUM_THREADS=$thread numactl -m 0 likwid-perfctr -m -g CLOCK -O -V 3 --stats -C 0-$((thread-1)) ../tradPower \
                -m ${matrixFolder}/${matrix} -c $thread -t 1 \
                -i $power -v $rcmFlag > tmp.txt
            cat tmp.txt >> ${raw_file}
            clock=$(grep "^Clock \[MHz\] STAT" tmp.txt | head -n ${selector} | tail -n 1 | cut -d"," -f 5)

            printf "%20s, %12d, %12d, %10.4f, %8.3f, %5d, %5d, %5d, %8.3f, %8.3f, %8.3f, %8.3f, %8.3f, %8.3f, %10.6f\n" ${matrix} ${nrows} ${nnz} ${perf} ${clock} ${thread} ${power} ${RCM} ${mem_data_trad} ${mem_bw_trad} ${l3_data_trad} ${l3_bw_trad} ${l2_data_trad} ${l2_bw_trad} ${preTime} >> ${tradPower_file}
        done
        let ctr=${ctr}+1
    done
fi


mtxPower_file="${folder}/mtxPower.csv"
printf "%20s, %12s, %12s, %10s, %10s, %8s, %8s, %5s, %5s, %5s, %5s, %5s, %5s, %8s, %8s, %8s, %8s, %8s, %8s, %8s, %8s, %8s, %8s, %8s, %8s, %10s, %10s\n" "# Matrix" "NROWS" "NNZ" "Perf" "TRAD_perf" "ClockRACE" "ClockTrad" "Thread" "Power" "CS" "N" "stage" "RCM" "MEM_race" "MEM_bw_race" "MEM_trad" "MEM_bw_trad" "L3_race" "L3_bw_race" "L3_trad" "L3_bw_trad" "L2_race" "L2_bw_race" "L2_trad" "L2_bw_trad" "RACE_Pre-time" "Total_Pre-time"> ${mtxPower_file}

powers=$(cut -d "," -f2 ${paramFileRACE} | tail -n +2)
cacheSizes=$(cut -d "," -f3 ${paramFileRACE} | tail -n +2)
RCMs=$(cut -d "," -f4 ${paramFileRACE} | tail -n +2)
stages=$(cut -d "," -f5 ${paramFileRACE} | tail -n +2)

ctr=1
for matrix in $matrix_name; do
    raw_file="${folder}/raw_mtxPower/${matrix}.txt"
    for thread in $threads; do
        power=$(echo ${powers} | cut -d" " -f${ctr})
        RCM_w_space=$(echo ${RCMs} | cut -d" " -f${ctr})
        RCM=$(echo ${RCM_w_space})
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

        RACE_CACHE_VIOLATION_CUTOFF=${safetyFactor} OMP_NUM_THREADS=$thread OMP_PROC_BIND=close OMP_PLACES=cores numactl -m 0 ../mtxPower \
            -m ${matrixFolder}/${matrix} -c $thread -t 1 \
            -i $power -n $node -s $cacheSizePerNode -v $rcmFlag &> tmp.txt

        cat tmp.txt >> ${raw_file}
        spmv_perf=$(readResult "tmp.txt" "SpMV power perf.")
        perf=$(readResult "tmp.txt" "RACE power perf.")
        nrows=$(cat tmp.txt | grep "Nrows =" | cut -d"=" -f2 | cut -d"," -f1)
        nnz=$(cat tmp.txt | grep "Nrows =" | cut -d"=" -f3 | cut -d"," -f1)
        racePreTime=$(cat tmp.txt | grep "RACE pre-processing time =" | cut -d"=" -f4 | cut -d"s" -f1)
        preTime=$(cat tmp.txt | grep "Total pre-processing time =" | cut -d"=" -f2 | cut -d"s" -f1)


        RACE_CACHE_VIOLATION_CUTOFF=${safetyFactor} OMP_NUM_THREADS=$thread numactl -m 0 likwid-perfctr -m -g MEM -O -V 3 --stats -C 0-$((thread-1)) ../mtxPower \
            -m ${matrixFolder}/${matrix} -c $thread -t 1 \
            -i $power -n $node -s $cacheSizePerNode -v $rcmFlag > tmp.txt
        cat tmp.txt >> ${raw_file}
        selector=1 #which region
        mem_data_trad=$(readDataVol "Memory data volume \[GBytes\] STAT" tmp.txt ${selector})
        mem_bw_trad=$(readMetric "Memory bandwidth \[MBytes/s\] STAT" tmp.txt ${selector})
        selector=2 #which region
        mem_data_race=$(readDataVol "Memory data volume \[GBytes\] STAT" tmp.txt ${selector})
        mem_bw_race=$(readMetric "Memory bandwidth \[MBytes/s\] STAT" tmp.txt ${selector})

        RACE_CACHE_VIOLATION_CUTOFF=${safetyFactor} OMP_NUM_THREADS=$thread numactl -m 0 likwid-perfctr -m -g L3 -O -V 3 --stats -C 0-$((thread-1)) ../mtxPower \
            -m ${matrixFolder}/${matrix} -c $thread -t 1 \
            -i $power -n $node -s $cacheSizePerNode -v $rcmFlag > tmp.txt
        cat tmp.txt >> ${raw_file}
        selector=1 #which region
        l3_data_trad=$(readDataVol "L3 data volume \[GBytes\] STAT" tmp.txt ${selector})
        l3_bw_trad=$(readMetric "L3 bandwidth \[MBytes/s\] STAT" tmp.txt ${selector})
        selector=2 #which region
        l3_data_race=$(readDataVol "L3 data volume \[GBytes\] STAT" tmp.txt ${selector})
        l3_bw_race=$(readMetric "L3 bandwidth \[MBytes/s\] STAT" tmp.txt ${selector})

        RACE_CACHE_VIOLATION_CUTOFF=${safetyFactor} OMP_NUM_THREADS=$thread numactl -m 0 likwid-perfctr -m -g L2 -O -V 3 --stats -C 0-$((thread-1)) ../mtxPower \
            -m ${matrixFolder}/${matrix} -c $thread -t 1 \
            -i $power -n $node -s $cacheSizePerNode -v $rcmFlag > tmp.txt
        cat tmp.txt >> ${raw_file}
        selector=1 #which region
        l2_data_trad=$(readDataVol "L2 data volume \[GBytes\] STAT" tmp.txt ${selector})
        l2_bw_trad=$(readMetric "L2 bandwidth \[MBytes/s\] STAT" tmp.txt ${selector})
        selector=2 #which region
        l2_data_race=$(readDataVol "L2 data volume \[GBytes\] STAT" tmp.txt ${selector})
        l2_bw_race=$(readMetric "L2 bandwidth \[MBytes/s\] STAT" tmp.txt ${selector})


        RACE_CACHE_VIOLATION_CUTOFF=${safetyFactor} OMP_NUM_THREADS=$thread numactl -m 0 likwid-perfctr -m -g CLOCK -O -V 3 --stats -C 0-$((thread-1)) ../mtxPower \
            -m ${matrixFolder}/${matrix} -c $thread -t 1 \
            -i $power -n $node -s $cacheSizePerNode -v $rcmFlag > tmp.txt
        cat tmp.txt >> ${raw_file}
        selector=1
        clock_trad=$(grep "^Clock \[MHz\] STAT" tmp.txt | head -n ${selector} | tail -n 1 | cut -d"," -f 5)
        selector=2
        clock_race=$(grep "^Clock \[MHz\] STAT" tmp.txt | head -n ${selector} | tail -n 1 | cut -d"," -f 5)

        #safetyFactor_wo_comma=$(echo ${safetyFactor} | sed -e "s@,@;@g")
        printf "%20s, %12d, %12d, %10.4f, %10.4f, %8.3f, %8.3f, %5d, %5d, %5.2f, %5d, %5d, %5d, %8.3f, %8.3f, %8.3f, %8.3f, %8.3f, %8.3f, %8.3f, %8.3f, %8.3f, %8.3f, %8.3f, %8.3f, %10.6f, %10.6f\n" ${matrix} ${nrows} ${nnz} ${perf} ${spmv_perf} ${clock_race} ${clock_trad} ${thread} ${power} ${cacheSizePerNode} ${node} ${stage} ${RCM} ${mem_data_race} ${mem_bw_race} ${mem_data_trad} ${mem_bw_trad} ${l3_data_race} ${l3_bw_race} ${l3_data_trad} ${l3_bw_trad} ${l2_data_race} ${l2_bw_race} ${l2_data_trad} ${l2_bw_trad} ${racePreTime} ${preTime} >> ${mtxPower_file}
    done
    let ctr=${ctr}+1
done

rm -rf tmp.txt
