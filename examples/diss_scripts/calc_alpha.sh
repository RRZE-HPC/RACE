#this script should benchmark coloring (RACE, ABMC, MC) and mklSymmSpMV executable

############# Settings ###################

#Read configurations from file
configFile=$1

matrixFolder=$(cat ${configFile} | grep "matrixFolder" | cut -d"=" -f2)
folder=$(cat ${configFile} | grep "folder" | cut -d"=" -f2)
threads=$(cat ${configFile} | grep "threads" | cut -d"=" -f2)
nodes=$(cat ${configFile} | grep "nodes" | cut -d"=" -f2)
RCMs=$(cat ${configFile} | grep "RCMs" | cut -d"=" -f2)
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

res_file="${folder}/spmv.csv"

rawFolder="${folder}/raw_spmv"
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

            KMP_WARNINGS=0 MKL_NUM_THREADS=$thread \
                OMP_NUM_THREADS=${threads} OMP_SCHEDULE=static \
                likwid-perfctr --stats -O -m -g MEM -C 0-$((threads-1)) ${execFolder}/spmv\
                -m "${matrixFolder}/${matrix}" -c ${thread} -t 1  -v\
                -T 1e-4 ${rcmFlag} > ${tmpFile}

            cat ${tmpFile} >> ${raw_file}

            if [[ $printHead == "1" ]]; then
                if [[ ${ctr} == 0 ]]; then
                    columns=$(printHeader ${tmpFile})
                    echo "Matrix,NROWS,NNZ,NNZR,NNZsymm,NNZRsymm,Thread,RCM,Iter,SpMV_perf,SpMV_datavol" > ${res_file}
                fi
            fi
            nrows_w_space=$(cat ${tmpFile} | grep "Nrows =" | cut -d"=" -f2 | cut -d"," -f1)
            nrows=$(echo ${nrows_w_space})
            nnz_w_space=$(cat ${tmpFile} | grep "NNZ =" | cut -d"=" -f3 | cut -d"," -f1)
            nnz=$(echo ${nnz_w_space})
            nnzr_w_space=$(cat ${tmpFile} | grep "NNZR =" | cut -d"=" -f4 | cut -d"," -f1)
            nnzr=$(echo ${nnzr_w_space})
            #for symmetric
            nnz_w_space=$(cat ${tmpFile} | grep "NNZ_symm =" | cut -d"=" -f5 | cut -d"," -f1)
            nnz_symm=$(echo ${nnz_w_space})
            nnzr_w_space=$(cat ${tmpFile} | grep "NNZR_symm =" | cut -d"=" -f6)
            nnzr_symm=$(echo ${nnzr_w_space})

            niter_w_space=$(cat ${tmpFile} | grep "Num iterations =" | cut -d"=" -f2)
            niter=$(echo ${niter_w_space})
            perf_w_space=$(cat ${tmpFile} | grep "PLAIN_SPMV :" | cut -d":" -f2 | cut -d"G" -f1)
            perf=$(echo ${perf_w_space})
            memVol_w_space=$(cat ${tmpFile} | grep "Memory data volume \[GBytes\] STAT" | cut -d"," -f2)
            memVol=$(echo ${memVol_w_space})
            echo "${matrix},${nrows},${nnz},${nnzr},${nnz_symm},${nnzr_symm},${thread},${RCM},${niter},${perf},${memVol}" >> ${res_file}

            let ctr=${ctr}+1
            rm -rf ${tmpFile}
        done
    done
done

#store aligned CSV file
sed 's/,/:,/g' ${res_file} | column -t -s: | sed 's/ ,/,/g' > "${alignFolder}/spmv.csv"

#compress raw folder
cd ${folder}
tar -cvzf raw_spmv.tar.gz raw_spmv
rm -rf raw_spmv
cd -

