#this script should benchmark coloring (RACE, ABMC, MC) and mklSymmSpMV executable

############# Settings ###################

#matrixFolder="/home/hk-project-benchfau/hd0705/matrices_diss"
#folder="symm_spmv_results"
#threads="38"
#nodes="1"
#RCMs="0 1"
#mkl="on"
#colorTypes="RACE ABMC MC"
#race_efficiencies="80,80 40"
#race_dists="1 2"
#execFolder="/home/hk-project-benchfau/hd0705/MatrixPower/examples/build_intel22"

#Read configurations from file
configFile=$1
matrixFolder=$(cat ${configFile} | grep "matrixFolder" | cut -d"=" -f2)
folder=$(cat ${configFile} | grep "folder" | cut -d"=" -f2)
threads=$(cat ${configFile} | grep "threads" | cut -d"=" -f2)
paramFileRACE=$(cat ${configFile} | grep "paramFileRACE" | cut -d"=" -f2)
paramFileTRAD=$(cat ${configFile} | grep "paramFileTRAD" | cut -d"=" -f2)
mkl=$(cat ${configFile} | grep "mkl" | cut -d"=" -f2)
colorTypes=$(cat ${configFile} | grep "colorTypes" | cut -d"=" -f2)
race_dists=2
execFolder=$(cat ${configFile} | grep "execFolder" | cut -d"=" -f2)
printHead=$(cat ${configFile} | grep "printHead" | cut -d"=" -f2)
############ Don't change anything below ############################

./check-state.sh config.txt

mkdir -p $folder
#align folder to store pretty CSV files
alignFolder="${folder}/align"
mkdir -p ${alignFolder}

matrix_name=$(cut -d "," -f1 ${paramFileRACE} | tail -n +2)

function readResult
{
    tmpFile=$1
    columns=$(grep "Obtained Perf of" ${tmpFile} | sed "s/Obtained Perf of//g" | cut -d":" -f2 | cut -d"G" -f 1)
    numColumns=$(echo $columns | grep -o "]" | wc -l)

    outStr=""
    for ((col_id=1; col_id<=${numColumns}; ++col_id)); do
        curCol=$(echo ${columns} | cut -d"]" -f${col_id} | sed 's/\[//')
        outStr="${outStr}${curCol},"
    done
    #remove last ',' and any spaces
    echo ${outStr} | sed 's/\(.*\),/\1 /' | sed 's/ //g'
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


function printHeader
{
    tmpFile=$1
    columns=$(grep "Obtained Perf of" ${tmpFile} | sed "s/Obtained Perf of//g" | cut -d":" -f1)
    outStr=""

    quartiles="0 25 50 75 100"
    #create header for each quartile
    for column in ${columns}; do
        for q in ${quartiles}; do
            outStr="${outStr}Perf_${column}_Q${q},"
        done
    done
    #remove last ,
    echo ${outStr} | sed 's/\(.*\),/\1 /'
}


if [[ ${mkl} == "on" ]]; then
    #TODO: adapt column number and so
    RCMs=$(cut -d "," -f2 ${paramFileTRAD} | tail -n +2)
    res_file="${folder}/mkl.csv"

    rawFolder="${folder}/raw_mkl"
    mkdir -p ${rawFolder}

    #get machine env
    ./machine-state.sh > ${rawFolder}/machine-state.txt


    ctr=0
    matCtr=1
    for matrix in $matrix_name; do
        raw_file="${rawFolder}/${matrix}.txt"
        for thread in $threads; do
            tmpFile="${rawFolder}/${matrix}.tmp"
            RCM_w_space=$(echo ${RCMs} | cut -d" " -f${matCtr})
            RCM=$(echo $RCM_w_space)
            rcmFlag=""
            if [[ $RCM == "1" ]]; then
                rcmFlag="-R"
            fi

            KMP_WARNINGS=0 MKL_NUM_THREADS=$thread OMP_NUM_THREADS=${thread} \
                OMP_SCHEDULE=static OMP_PROC_BIND=close OMP_PLACES=cores \
                taskset -c 0-$((thread-1)) ${execFolder}/mklSymmSpMV \
                -m "${matrixFolder}/${matrix}" -c ${thread} -T 1e-4 ${rcmFlag} > ${tmpFile}

            cat ${tmpFile} >> ${raw_file}
            preTime_w_space=$(cat ${tmpFile} | grep "Total pre-processing time =" | cut -d"=" -f2 | cut -d"s" -f1)
            preTime=$(echo ${preTime_w_space})
            perf=$(readResult "${tmpFile}") #TODO: select only SymmSpMV
            nrows=$(cat ${tmpFile} | grep "NROWS =" | cut -d"=" -f2 | cut -d"," -f1)
            nnz=$(cat ${tmpFile} | grep "NROWS =" | cut -d"=" -f3 | cut -d"," -f1)

            selector=3
            KMP_WARNINGS=0 MKL_NUM_THREADS=$thread OMP_NUM_THREADS=${thread} \
                OMP_SCHEDULE=static OMP_PROC_BIND=close OMP_PLACES=cores \
                taskset -c 0-$((thread-1)) likwid-perfctr -m -g MEM -O --stats -c 0-$((thread-1)) ${execFolder}/mklSymmSpMV \
                -m "${matrixFolder}/${matrix}" -c ${thread} -T 1e-4 ${rcmFlag} > ${tmpFile}

            cat ${tmpFile} >> ${raw_file}
            mem_data_trad=$(readDataVol "Memory data volume \[GBytes\] STAT" ${tmpFile} ${selector})
            mem_bw_trad=$(readMetric "Memory bandwidth \[MBytes/s\] STAT" ${tmpFile} ${selector})


            KMP_WARNINGS=0 MKL_NUM_THREADS=$thread OMP_NUM_THREADS=${thread} \
                OMP_SCHEDULE=static OMP_PROC_BIND=close OMP_PLACES=cores \
                taskset -c 0-$((thread-1)) likwid-perfctr -m -g L3 -O --stats -c 0-$((thread-1)) ${execFolder}/mklSymmSpMV \
                -m "${matrixFolder}/${matrix}" -c ${thread} -T 1e-4 ${rcmFlag} > ${tmpFile}

            cat ${tmpFile} >> ${raw_file}
            l3_data_trad=$(readDataVol "L3 data volume \[GBytes\] STAT" ${tmpFile} ${selector})
            l3_bw_trad=$(readMetric "L3 bandwidth \[MBytes/s\] STAT" ${tmpFile} ${selector})


            KMP_WARNINGS=0 MKL_NUM_THREADS=$thread OMP_NUM_THREADS=${thread} \
                OMP_SCHEDULE=static OMP_PROC_BIND=close OMP_PLACES=cores \
                taskset -c 0-$((thread-1)) likwid-perfctr -m -g L2 -O --stats -c 0-$((thread-1)) ${execFolder}/mklSymmSpMV \
                -m "${matrixFolder}/${matrix}" -c ${thread} -T 1e-4 ${rcmFlag} > ${tmpFile}

            cat ${tmpFile} >> ${raw_file}
            l2_data_trad=$(readDataVol "L2 data volume \[GBytes\] STAT" ${tmpFile} ${selector})
            l2_bw_trad=$(readMetric "L2 bandwidth \[MBytes/s\] STAT" ${tmpFile} ${selector})

            KMP_WARNINGS=0 MKL_NUM_THREADS=$thread OMP_NUM_THREADS=${thread} \
                OMP_SCHEDULE=static OMP_PROC_BIND=close OMP_PLACES=cores \
                taskset -c 0-$((thread-1)) likwid-perfctr -m -g CLOCK -O --stats -c 0-$((thread-1)) ${execFolder}/mklSymmSpMV \
                -m "${matrixFolder}/${matrix}" -c ${thread} -T 1e-4 ${rcmFlag} > ${tmpFile}

            cat ${tmpFile} >> ${raw_file}
            clock=$(grep "^Clock \[MHz\] STAT" ${tmpFile} | head -n ${selector} | tail -n 1 | cut -d"," -f 5)

            if [[ $printHead == "1" ]]; then
                if [[ ${ctr} == 0 ]]; then
                    columns=$(printHeader ${tmpFile})
                    printf "%20s, %12s, %12s, %8s, %5s, %5s, %8s, %8s, %8s, %8s, %8s, %8s, %10s, %s\n" "#Matrix" "NROWS" "NNZ" "Clock" "Thread" "RCM" "MEM" "MEM_bw" "L3" "L3_bw" "L2" "L2_bw" "Total_Pre-time" ${columns} > ${res_file}
                fi
            fi

            printf "%20s, %12d, %12d, %8.3f, %5d, %5d, %8.3f, %8.3f, %8.3f, %8.3f, %8.3f, %8.3f, %10.6f, %s\n" ${matrix} ${nrows} ${nnz} ${clock} ${thread} ${RCM} ${mem_data_trad} ${mem_bw_trad} ${l3_data_trad} ${l3_bw_trad} ${l2_data_trad} ${l2_bw_trad} ${preTime} "${perf}" >> ${res_file}

            let ctr=${ctr}+1
            rm -rf ${tmpFile}
        done
        let matCtr=${matCtr}+1
    done

    #store aligned CSV file
    sed 's/,/:,/g' ${res_file} | column -t -s: | sed 's/ ,/,/g' > "${alignFolder}/mkl.csv"

    #compress raw folder
    cd ${folder}
    tar -cvzf raw_mkl.tar.gz raw_mkl
    rm -rf raw_mkl
    cd -
fi

#TODO: adapt column number and so
RCMs=$(cut -d "," -f2 ${paramFileRACE} | tail -n +2)
race_efficiencies=$(cut -d "," -f3 ${paramFileRACE} | tail -n +2)


for colorType in ${colorTypes}; do
    res_file="${folder}/${colorType}.csv"

    rawFolder="${folder}/raw_${colorType}"
    mkdir -p ${rawFolder}
    #get machine env
    ./machine-state.sh > ${rawFolder}/machine-state.txt

    ctr=0
    matCtr=1
    for matrix in $matrix_name; do
        raw_file="${rawFolder}/${matrix}.txt"
        for thread in $threads; do
            RCM_w_space=$(echo ${RCMs} | cut -d" " -f${matCtr})
            RCM=$(echo $RCM_w_space)

            eff=40
            if [[ ${colorType} == "RACE" ]]; then #pinning left to RACE
                eff_w_space=$(echo ${race_efficiencies} | cut -d" " -f${matCtr})
                eff=$(echo ${eff_w_space} | sed 's/+/,/g')
            fi

            tmpFile="${rawFolder}/${matrix}.tmp"
            rcmFlag=""
            if [[ $RCM == "1" ]]; then
                rcmFlag="-R"
            fi

            execArgs="-m ${matrixFolder}/${matrix} -c ${thread} -t 1  -v -p FILL -T 1e-4 -C ${colorType} ${rcmFlag}"
            if [[ ${colorType} == "RACE" ]]; then #pinning left to RACE
                KMP_WARNINGS=0 MKL_NUM_THREADS=$thread \
                    OMP_NUM_THREADS=${thread} OMP_SCHEDULE=static \
                    COLOR_DISTANCE=${dist} RACE_EFFICIENCY=${eff} \
                    taskset -c 0-$((thread-1)) ${execFolder}/colorSymmSpMV ${execArgs} > ${tmpFile}
            else #here pinning via OMP
                 KMP_WARNINGS=0 MKL_NUM_THREADS=$thread \
                    OMP_NUM_THREADS=${threads} OMP_SCHEDULE=static \
                    OMP_PROC_BIND=close OMP_PLACES=cores \
                    COLOR_DISTANCE=${dist} RACE_EFFICIENCY=${eff} \
                    taskset -c 0-$((thread-1)) ${execFolder}/colorSymmSpMV ${execArgs} > ${tmpFile}
            fi
            cat ${tmpFile} >> ${raw_file}
            niter_w_space=$(cat ${tmpFile} | grep "Num iterations =" | cut -d"=" -f2)
            niter=$(echo ${niter_w_space})
            preTime_w_space=$(cat ${tmpFile} | grep "Total pre-processing time =" | cut -d"=" -f2 | cut -d"s" -f1)
            preTime=$(echo ${preTime_w_space})

            racePreTime_w_space=$(cat ${tmpFile} | grep "RACE pre-processing time =" | cut -d"=" -f2 | cut -d"s" -f1)
            racePreTime=$(echo ${racePreTime_w_space})


            perf=$(readResult "${tmpFile}")
            nrows=$(cat ${tmpFile} | grep "NROWS =" | cut -d"=" -f2 | cut -d"," -f1)
            nnz=$(cat ${tmpFile} | grep "NROWS =" | cut -d"=" -f3 | cut -d"," -f1)

            if [[ ${colorType} == "RACE" ]]; then #pinning left to RACE
                KMP_WARNINGS=0 MKL_NUM_THREADS=$thread \
                    OMP_NUM_THREADS=${thread} OMP_SCHEDULE=static \
                    COLOR_DISTANCE=${dist} RACE_EFFICIENCY=${eff} \
                    taskset -c 0-$((thread-1)) likwid-perfctr -m -g MEM -O --stats -c 0-$((thread-1)) ${execFolder}/colorSymmSpMV ${execArgs} > ${tmpFile}
            else #here pinning via OMP
                 KMP_WARNINGS=0 MKL_NUM_THREADS=$thread \
                    OMP_NUM_THREADS=${threads} OMP_SCHEDULE=static \
                    OMP_PROC_BIND=close OMP_PLACES=cores \
                    COLOR_DISTANCE=${dist} RACE_EFFICIENCY=${eff} \
                    taskset -c 0-$((thread-1)) likwid-perfctr -m -g MEM -O --stats -c 0-$((thread-1)) ${execFolder}/colorSymmSpMV ${execArgs} > ${tmpFile}
            fi
            cat ${tmpFile} >> ${raw_file}
            selector=1
            mem_data_trad=$(readDataVol "Memory data volume \[GBytes\] STAT" ${tmpFile} ${selector})
            mem_bw_trad=$(readMetric "Memory bandwidth \[MBytes/s\] STAT" ${tmpFile} ${selector})
            selector=3
            mem_data_race=$(readDataVol "Memory data volume \[GBytes\] STAT" ${tmpFile} ${selector})
            mem_bw_race=$(readMetric "Memory bandwidth \[MBytes/s\] STAT" ${tmpFile} ${selector})

            if [[ ${colorType} == "RACE" ]]; then #pinning left to RACE
                KMP_WARNINGS=0 MKL_NUM_THREADS=$thread \
                    OMP_NUM_THREADS=${thread} OMP_SCHEDULE=static \
                    COLOR_DISTANCE=${dist} RACE_EFFICIENCY=${eff} \
                    taskset -c 0-$((thread-1)) likwid-perfctr -m -g L3 -O --stats -c 0-$((thread-1)) ${execFolder}/colorSymmSpMV ${execArgs} > ${tmpFile}
            else #here pinning via OMP
                 KMP_WARNINGS=0 MKL_NUM_THREADS=$thread \
                    OMP_NUM_THREADS=${threads} OMP_SCHEDULE=static \
                    OMP_PROC_BIND=close OMP_PLACES=cores \
                    COLOR_DISTANCE=${dist} RACE_EFFICIENCY=${eff} \
                    taskset -c 0-$((thread-1)) likwid-perfctr -m -g L3 -O --stats -c 0-$((thread-1)) ${execFolder}/colorSymmSpMV ${execArgs} > ${tmpFile}
            fi
            cat ${tmpFile} >> ${raw_file}
            selector=1
            l3_data_trad=$(readDataVol "L3 data volume \[GBytes\] STAT" ${tmpFile} ${selector})
            l3_bw_trad=$(readMetric "L3 bandwidth \[MBytes/s\] STAT" ${tmpFile} ${selector})
            selector=3
            l3_data_race=$(readDataVol "L3 data volume \[GBytes\] STAT" ${tmpFile} ${selector})
            l3_bw_race=$(readMetric "L3 bandwidth \[MBytes/s\] STAT" ${tmpFile} ${selector})

            if [[ ${colorType} == "RACE" ]]; then #pinning left to RACE
                KMP_WARNINGS=0 MKL_NUM_THREADS=$thread \
                    OMP_NUM_THREADS=${thread} OMP_SCHEDULE=static \
                    COLOR_DISTANCE=${dist} RACE_EFFICIENCY=${eff} \
                    taskset -c 0-$((thread-1)) likwid-perfctr -m -g L2 -O --stats -c 0-$((thread-1)) ${execFolder}/colorSymmSpMV ${execArgs} > ${tmpFile}
            else #here pinning via OMP
                 KMP_WARNINGS=0 MKL_NUM_THREADS=$thread \
                    OMP_NUM_THREADS=${threads} OMP_SCHEDULE=static \
                    OMP_PROC_BIND=close OMP_PLACES=cores \
                    COLOR_DISTANCE=${dist} RACE_EFFICIENCY=${eff} \
                    taskset -c 0-$((thread-1)) likwid-perfctr -m -g L2 -O --stats -c 0-$((thread-1)) ${execFolder}/colorSymmSpMV ${execArgs} > ${tmpFile}
            fi
            cat ${tmpFile} >> ${raw_file}
            selector=1
            l2_data_trad=$(readDataVol "L2 data volume \[GBytes\] STAT" ${tmpFile} ${selector})
            l2_bw_trad=$(readMetric "L2 bandwidth \[MBytes/s\] STAT" ${tmpFile} ${selector})
            selector=3
            l2_data_race=$(readDataVol "L2 data volume \[GBytes\] STAT" ${tmpFile} ${selector})
            l2_bw_race=$(readMetric "L2 bandwidth \[MBytes/s\] STAT" ${tmpFile} ${selector})

            if [[ ${colorType} == "RACE" ]]; then #pinning left to RACE
                KMP_WARNINGS=0 MKL_NUM_THREADS=$thread \
                    OMP_NUM_THREADS=${thread} OMP_SCHEDULE=static \
                    COLOR_DISTANCE=${dist} RACE_EFFICIENCY=${eff} \
                    taskset -c 0-$((thread-1)) likwid-perfctr -m -g CLOCK -O --stats -c 0-$((thread-1)) ${execFolder}/colorSymmSpMV ${execArgs} > ${tmpFile}
            else #here pinning via OMP
                 KMP_WARNINGS=0 MKL_NUM_THREADS=$thread \
                    OMP_NUM_THREADS=${threads} OMP_SCHEDULE=static \
                    OMP_PROC_BIND=close OMP_PLACES=cores \
                    COLOR_DISTANCE=${dist} RACE_EFFICIENCY=${eff} \
                    taskset -c 0-$((thread-1)) likwid-perfctr -m -g CLOCK -O --stats -c 0-$((thread-1)) ${execFolder}/colorSymmSpMV ${execArgs} > ${tmpFile}
            fi
            cat ${tmpFile} >> ${raw_file}
            selector=1
            clock_trad=$(grep "^Clock \[MHz\] STAT" ${tmpFile} | head -n ${selector} | tail -n 1 | cut -d"," -f 5)
            selector=3
            clock_race=$(grep "^Clock \[MHz\] STAT" ${tmpFile} | head -n ${selector} | tail -n 1 | cut -d"," -f 5)

            if [[ $printHead == "1" ]]; then
                if [[ ${ctr} == 0 ]]; then
                    columns=$(printHeader ${tmpFile})
                    if [[ ${colorType} == "RACE" ]]; then
                        printf "%20s, %12s, %12s, %12s, %8s, %8s, %5s, %5s, %8s, %8s, %8s, %8s, %8s, %8s, %8s, %8s, %8s, %8s, %8s, %8s, %10s, %10s, %s\n" "# Matrix" "NROWS" "NNZ" "RACE_eff" "ClockRACE" "ClockTrad" "Thread" "RCM" "MEM_race" "MEM_bw_race" "MEM_trad" "MEM_bw_trad" "L3_race" "L3_bw_race" "L3_trad" "L3_bw_trad" "L2_race" "L2_bw_race" "L2_trad" "L2_bw_trad" "RACE_Pre-time" "Total_Pre-time" ${columns}> ${res_file}
                    else
                        printf "%20s, %12s, %12s, %8s, %8s, %5s, %5s, %8s, %8s, %8s, %8s, %8s, %8s, %8s, %8s, %8s, %8s, %8s, %8s, %10s, %10s, %s\n" "# Matrix" "NROWS" "NNZ" "ClockRACE" "ClockTrad" "Thread" "RCM" "MEM_race" "MEM_bw_race" "MEM_trad" "MEM_bw_trad" "L3_race" "L3_bw_race" "L3_trad" "L3_bw_trad" "L2_race" "L2_bw_race" "L2_trad" "L2_bw_trad" "RACE_Pre-time" "Total_Pre-time" ${columns}> ${res_file}
                    fi
                fi
            fi
            if [[ ${colorType} == "RACE" ]]; then
                effPrintStr=$(echo ${eff} | sed 's/,/+/g')
                printf "%20s, %12d, %12d, %12s, %8.3f, %8.3f, %5d, %5d, %8.3f, %8.3f, %8.3f, %8.3f, %8.3f, %8.3f, %8.3f, %8.3f, %8.3f, %8.3f, %8.3f, %8.3f, %10.6f, %10.6f, %s\n" ${matrix} ${nrows} ${nnz} ${effPrintStr} ${clock_race} ${clock_trad} ${thread} ${RCM} ${mem_data_race} ${mem_bw_race} ${mem_data_trad} ${mem_bw_trad} ${l3_data_race} ${l3_bw_race} ${l3_data_trad} ${l3_bw_trad} ${l2_data_race} ${l2_bw_race} ${l2_data_trad} ${l2_bw_trad} ${racePreTime} ${preTime} "${perf}" >> ${res_file}
            else
                printf "%20s, %12d, %12d, %8.3f, %8.3f, %5d, %5d, %8.3f, %8.3f, %8.3f, %8.3f, %8.3f, %8.3f, %8.3f, %8.3f, %8.3f, %8.3f, %8.3f, %8.3f, %10.6f, %10.6f, %s\n" ${matrix} ${nrows} ${nnz} ${clock_race} ${clock_trad} ${thread} ${RCM} ${mem_data_race} ${mem_bw_race} ${mem_data_trad} ${mem_bw_trad} ${l3_data_race} ${l3_bw_race} ${l3_data_trad} ${l3_bw_trad} ${l2_data_race} ${l2_bw_race} ${l2_data_trad} ${l2_bw_trad} ${racePreTime} ${preTime} "${perf}" >> ${res_file}
            fi

            let ctr=${ctr}+1
            rm -rf ${tmpFile}
        done
        let matCtr=${matCtr}+1
    done

    #store aligned CSV file
    sed 's/,/:,/g' ${res_file} | column -t -s: | sed 's/ ,/,/g' > "${alignFolder}/${colorType}.csv"

    #compress raw folder
    cd ${folder}
    tar -cvzf raw_${colorType}.tar.gz raw_${colorType}
    rm -rf raw_${colorType}
    cd -

done
