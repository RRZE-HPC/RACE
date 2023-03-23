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
#execFolder="/home/hk-project-benchfau/hd0705/MatrixPower/examples/build_intel22"

#Read configurations from file
configFile=$1

matrixFolder=$(cat ${configFile} | grep "matrixFolder" | cut -d"=" -f2)
folder=$(cat ${configFile} | grep "folder" | cut -d"=" -f2)
threads=$(cat ${configFile} | grep "threads" | cut -d"=" -f2)
nodes=$(cat ${configFile} | grep "nodes" | cut -d"=" -f2)
mkl=$(cat ${configFile} | grep "mkl" | cut -d"=" -f2)
execFolder=$(cat ${configFile} | grep "execFolder" | cut -d"=" -f2)
printHead=$(cat ${configFile} | grep "printHead" | cut -d"=" -f2)
colorType=$(cat ${configFile} | grep "colorType" | cut -d"=" -f2) #only RACE supported
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

res_file="${folder}/${colorType}.csv"
convFolder="${folder}/convergence/${colorType}"
mkdir -p ${convFolder}
rawFolder="${folder}/raw_${colorType}"
mkdir -p ${rawFolder}
#get machine env
./machine-state.sh > ${rawFolder}/machine-state.txt

inCtr=0

PARAM_file="${folder}/kacz_best_parameters/${colorType}.csv"
echo $PARAM_file
fileLen=$(wc -l ${PARAM_file} | cut -d" " -f1)
numMatrices=$(echo "${fileLen}-1" | bc -l)
for(( ctr=0; ctr < ${numMatrices}; ctr=${ctr}+1 )); do
    #read RCM and efficiency from file
    line=$(head -n $((ctr+2)) ${PARAM_file} | tail -n 1)
    PARAM_matrix_w_space=$(echo ${line} | cut -d"," -f2)
    PARAM_matrix=$(echo ${PARAM_matrix_w_space})
    matrix=${PARAM_matrix}
    tmpFile="${rawFolder}/${matrix}.tmp"
    raw_file="${rawFolder}/${matrix}.txt"
    PARAMrcm_w_space=$(echo ${line} | cut -d"," -f3)
    PARAMrcm=$(echo ${PARAMrcm_w_space})

    PARAMiter_w_space=$(echo ${line} | cut -d"," -f4) #this is input iteration
    PARAMiter=$(echo ${PARAMiter_w_space})
    iter=$(echo "${PARAMiter}" | bc -l) #scale by 10, to allow incase of 10x bad convergence

    PARAMerr_w_space=$(echo ${line} | cut -d"," -f5)
    PARAMerr=$(echo ${PARAMerr_w_space})

    eff="40"
    if [[ ${colorType} == "RACE" ]]; then
        PARAMeff_w_space==$(echo ${line} | cut -d"," -f6)
        PARAMeff=$(echo ${PARAMeff_w_space})
        PARAMeff_replaced=$(echo ${PARAMeff} | sed 's/+/,/g')
        eff=${PARAMeff_replaced}
    fi
    #not differentiating for any coloring methods the parameters,
    #because we compare with RACE the other methods and
    #reordering in RACE can change the errNorm and/or iter
    echo $PARAMrcm
    rcmFlag=""
    if [[ $PARAMrcm == "1" ]]; then
        rcmFlag="-R"
    fi

    for thread in $threads; do
        if [[ ${colorType} == "RACE" ]]; then
            #pinning left to RACE
            #iterations automatically decide
            echo "KMP_WARNINGS=0 MKL_NUM_THREADS=$thread OMP_NUM_THREADS=${thread} OMP_SCHEDULE=static COLOR_DISTANCE=2 RACE_EFFICIENCY=${eff} taskset -c 0-$((thread-1)) ${execFolder}/colorKACZ -m ${matrixFolder}/${matrix} -c ${thread} -t 1  -p FILL -C ${colorType} ${rcmFlag} -i ${iter} -e ${PARAMerr} -v -f ${convFolder}/${matrix}.txt"
            KMP_WARNINGS=0 MKL_NUM_THREADS=$thread \
                OMP_NUM_THREADS=${thread} OMP_SCHEDULE=static \
                COLOR_DISTANCE=2 RACE_EFFICIENCY=${eff} \
                taskset -c 0-$((thread-1)) ${execFolder}/colorKACZ \
                -m "${matrixFolder}/${matrix}" -c ${thread} -t 1  -p FILL \
                -C ${colorType} ${rcmFlag} \
                -i ${iter} -e ${PARAMerr} \
                -v -f "${convFolder}/${matrix}_thread_${thread}.txt" > ${tmpFile}
        elif [[ ${colorType} == "SERIAL" ]]; then
            #here pinning via OMP
            #try to achieve the same error as RACE
            #and give 2*iterations of RACE as max. iter
            KMP_WARNINGS=0 MKL_NUM_THREADS=1 \
                OMP_NUM_THREADS=1 OMP_SCHEDULE=static \
                OMP_PROC_BIND=close OMP_PLACES=cores \
                COLOR_DISTANCE=2 RACE_EFFICIENCY=${eff} \
                taskset -c 0-$((thread-1)) ${execFolder}/serialKACZ \
                -m "${matrixFolder}/${matrix}" -c ${thread} -t 1  -p FILL \
                -C ${colorType} ${rcmFlag} \
                -i ${iter} -e ${PARAMerr} \
                -v -f "${convFolder}/${matrix}_thread_${thread}.txt" > ${tmpFile}
        else
            #here pinning via OMP
            #try to achieve the same error as RACE
            #and give 2*iterations of RACE as max. iter
            KMP_WARNINGS=0 MKL_NUM_THREADS=$thread \
                OMP_NUM_THREADS=${thread} OMP_SCHEDULE=static \
                OMP_PROC_BIND=close OMP_PLACES=cores \
                COLOR_DISTANCE=2 RACE_EFFICIENCY=${eff} \
                taskset -c 0-$((thread-1)) ${execFolder}/colorKACZ \
                -m "${matrixFolder}/${matrix}" -c ${thread} -t 1  -p FILL \
                -C ${colorType} ${rcmFlag} \
                -i ${iter} -e ${PARAMerr} \
                -v -f "${convFolder}/${matrix}_thread_${thread}.txt" > ${tmpFile}
        fi

        echo "Dumping raw file for Matrix=${matrix}, thread=${thread}" >> ${raw_file}
        cat ${tmpFile} >> ${raw_file}

        if [[ $printHead == "1" ]]; then
            if [[ ${inCtr} == 0 ]]; then
                columns=$(printHeader ${tmpFile})
                if [[ ${colorType} == "RACE" ]]; then
                    echo "Matrix,Thread,RCM,RACE_efficiency,Iter,ResNorm,ErrNorm,ActualIter,Preprocessing-time,${columns}" > ${res_file}
                else
                    echo "Matrix,Thread,RCM,Iter,ResNorm,ErrNorm,ActualIter,Preprocessing-time,${columns}" > ${res_file}
                fi
            fi
        fi
        niter_w_space=$(cat ${tmpFile} | grep "Num iterations =" | cut -d"=" -f2)
        niter=$(echo ${niter_w_space})
        preTime_w_space=$(cat ${tmpFile} | grep "Total pre-processing time =" | cut -d"=" -f2 | cut -d"s" -f1)
        preTime=$(echo ${preTime_w_space})
        perfRes=$(readResult "${tmpFile}")
        effPrintStr=$(echo ${eff} | sed 's/,/+/g')
        resNorm=$(grep "Convergence results:" ${tmpFile} | cut -d"=" -f2 | cut -d"," -f1)
        errNorm=$(grep "Convergence results:" ${tmpFile} | cut -d"=" -f3 | cut -d"," -f1)
        actualIter=$(grep "Convergence results:" ${tmpFile} | cut -d"=" -f4 | cut -d"," -f1)

        if [[ ${colorType} == "RACE" ]]; then
            echo "${matrix},${thread},${RCM},${effPrintStr},${niter},${resNorm},${errNorm},${actualIter},${preTime},${perfRes}" >> ${res_file}
        else
            echo "${matrix},${thread},${RCM},${niter},${resNorm},${errNorm},${actualIter},${preTime},${perfRes}" >> ${res_file}
        fi

        rm -rf ${tmpFile}

        let inCtr=${inCtr}+1
        #done
    done
done

#store aligned CSV file
sed 's/,/:,/g' ${res_file} | column -t -s: | sed 's/ ,/,/g' > "${alignFolder}/${colorType}.csv"

#compress raw folder
cd ${folder}
tar -cvzf raw_${colorType}.tar.gz raw_${colorType}
rm -rf raw_${colorType}
cd -

