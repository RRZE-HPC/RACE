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
nodes=$(cat ${configFile} | grep "nodes" | cut -d"=" -f2)
RCMs=$(cat ${configFile} | grep "RCMs" | cut -d"=" -f2)
mkl=$(cat ${configFile} | grep "mkl" | cut -d"=" -f2)
colorTypes=$(cat ${configFile} | grep "colorTypes" | cut -d"=" -f2)
race_efficiencies=$(cat ${configFile} | grep "race_efficiencies" | cut -d"=" -f2)
race_dists=$(cat ${configFile} | grep "race_dists" | cut -d"=" -f2)
execFolder=$(cat ${configFile} | grep "execFolder" | cut -d"=" -f2)
printHead=$(cat ${configFile} | grep "printHead" | cut -d"=" -f2)
############ Don't change anything below ############################

./check-state.sh config.txt

mkdir -p $folder
#align folder to store pretty CSV files
alignFolder="${folder}/align"
mkdir -p ${alignFolder}

#get matrix names
cd $matrixFolder
matrix_name=

while read -r i
do
    matrix_name=$matrix_name" "$i
done < <(find *.mtx)

echo $matrix_name

cd -


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
    res_file="${folder}/mkl.csv"

    rawFolder="${folder}/raw_mkl"
    mkdir -p ${rawFolder}

    #get machine env
    ./machine-state.sh > ${rawFolder}/machine-state.txt


    ctr=0
    for matrix in $matrix_name; do
        raw_file="${rawFolder}/${matrix}.txt"
        for thread in $threads; do
            for RCM in $RCMs; do
                tmpFile="${rawFolder}/${matrix}.tmp"
                rcmFlag=""
                if [[ $RCM == "1" ]]; then
                    rcmFlag="-R"
                fi

                KMP_WARNINGS=0 MKL_NUM_THREADS=$thread OMP_NUM_THREADS=${thread} \
                OMP_SCHEDULE=static OMP_PROC_BIND=close OMP_PLACES=cores \
                taskset -c 0-$((thread-1)) ${execFolder}/mklSymmSpMV \
                -m "${matrixFolder}/${matrix}" -c ${thread} -T 1e-4 ${rcmFlag} > ${tmpFile}

                cat ${tmpFile} >> ${raw_file}
                if [[ $printHead == "1" ]]; then
                    if [[ ${ctr} == 0 ]]; then
                        columns=$(printHeader ${tmpFile})
                        echo "Matrix,Thread,RCM,Iter,Preprocessing-time,${columns}" > ${res_file}
                    fi
                fi
                niter_w_space=$(cat ${tmpFile} | grep "Num iterations =" | cut -d"=" -f2)
                niter=$(echo ${niter_w_space})
                preTime_w_space=$(cat ${tmpFile} | grep "Total pre-processing time =" | cut -d"=" -f2 | cut -d"s" -f1)
                preTime=$(echo ${preTime_w_space})
                perfRes=$(readResult "${tmpFile}")
                echo "${matrix},${thread},${RCM},${niter},${preTime},${perfRes}" >> ${res_file}
                let ctr=${ctr}+1
                rm -rf ${tmpFile}
            done
        done
    done

    #store aligned CSV file
    sed 's/,/:,/g' ${res_file} | column -t -s: | sed 's/ ,/,/g' > "${alignFolder}/mkl.csv"

    #compress raw folder
    cd ${folder}
    tar -cvzf raw_mkl.tar.gz raw_mkl
    rm -rf raw_mkl
    cd -
fi


for colorType in ${colorTypes}; do
    res_file="${folder}/${colorType}.csv"

    rawFolder="${folder}/raw_${colorType}"
    mkdir -p ${rawFolder}
    #get machine env
    ./machine-state.sh > ${rawFolder}/machine-state.txt

    ctr=0
    for matrix in $matrix_name; do
        raw_file="${rawFolder}/${matrix}.txt"
        for thread in $threads; do
            for RCM in $RCMs; do
                race_effs_parse="40"
                dists_parse="2"
                if [[ ${colorType} == "RACE" ]]; then
                    race_effs_parse=${race_efficiencies}
                    dists_parse=${race_dists}
                fi
                for eff in ${race_effs_parse}; do
                    for dist in ${dists_parse}; do
                        tmpFile="${rawFolder}/${matrix}.tmp"
                        rcmFlag=""
                        if [[ $RCM == "1" ]]; then
                            rcmFlag="-R"
                        fi

                        if [[ ${colorType} == "RACE" ]]; then #pinning left to RACE
                            KMP_WARNINGS=0 MKL_NUM_THREADS=$thread \
                                OMP_NUM_THREADS=${thread} OMP_SCHEDULE=static \
                                COLOR_DISTANCE=${dist} RACE_EFFICIENCY=${eff} \
                                taskset -c 0-$((thread-1)) ${execFolder}/colorSymmSpMV \
                                -m "${matrixFolder}/${matrix}" -c ${thread} -t 1  -v -p FILL \
                                -T 1e-4 -C ${colorType} ${rcmFlag} > ${tmpFile}
                        else #here pinning via OMP
                            KMP_WARNINGS=0 MKL_NUM_THREADS=$thread \
                                OMP_NUM_THREADS=${threads} OMP_SCHEDULE=static \
                                OMP_PROC_BIND=close OMP_PLACES=cores \
                                COLOR_DISTANCE=${dist} RACE_EFFICIENCY=${eff} \
                                taskset -c 0-$((thread-1)) ${execFolder}/colorSymmSpMV \
                                -m "${matrixFolder}/${matrix}" -c ${thread} -t 1  -v -p FILL \
                                -T 1e-4 -C ${colorType} ${rcmFlag} > ${tmpFile}
                        fi

                        cat ${tmpFile} >> ${raw_file}

                        if [[ $printHead == "1" ]]; then
                            if [[ ${ctr} == 0 ]]; then
                                columns=$(printHeader ${tmpFile})
                                if [[ ${colorType} == "RACE" ]]; then
                                    echo "Matrix,Thread,RCM,RACE_efficiency,RACE_distance,Iter,Preprocessing-time,${columns}" > ${res_file}
                                else
                                    echo "Matrix,Thread,RCM,Iter,Preprocessing-time,${columns}" > ${res_file}
                                fi
                            fi
                        fi
                        niter_w_space=$(cat ${tmpFile} | grep "Num iterations =" | cut -d"=" -f2)
                        niter=$(echo ${niter_w_space})
                        preTime_w_space=$(cat ${tmpFile} | grep "Total pre-processing time =" | cut -d"=" -f2 | cut -d"s" -f1)
                        preTime=$(echo ${preTime_w_space})
                        perfRes=$(readResult "${tmpFile}")
                        if [[ ${colorType} == "RACE" ]]; then
                            effPrintStr=$(echo ${eff} | sed 's/,/+/g')
                            echo "${matrix},${thread},${RCM},${effPrintStr},${dist},${niter},${preTime},${perfRes}" >> ${res_file}
                        else
                            echo "${matrix},${thread},${RCM},${niter},${preTime},${perfRes}" >> ${res_file}
                        fi

                        let ctr=${ctr}+1
                        rm -rf ${tmpFile}
                    done
                done
            done
        done
    done

    #store aligned CSV file
    sed 's/,/:,/g' ${res_file} | column -t -s: | sed 's/ ,/,/g' > "${alignFolder}/${colorType}.csv"

    #compress raw folder
    cd ${folder}
    tar -cvzf raw_${colorType}.tar.gz raw_${colorType}
    rm -rf raw_${colorType}
    cd -

done
