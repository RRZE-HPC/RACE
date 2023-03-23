#this script should benchmark coloring (RACE, ABMC, MC) and mklSymmSpMV executable

############# Settings ###################

#matrixFolder="/home/hk-project-benchfau/hd0705/matrices_diss"
#folder="symm_spmv_results"
#threads="38"
#nodes="1"
#mkl="on"
#colorTypes="RACE ABMC MC"
#race_efficiencies="80,80 40"
#race_dists="1 2"
#execFolder="/home/hk-project-benchfau/hd0705/MatrixPower/examples/build_intel22"

#Read configurations from file
configFile=$1

matrixFolder=$(cat ${configFile} | grep "matrixFolder" | cut -d"=" -f2)
folder=$(cat ${configFile} | grep "folder" | cut -d"=" -f2)
thread=$(cat ${configFile} | grep "thread" | cut -d"=" -f2)
nodes=$(cat ${configFile} | grep "nodes" | cut -d"=" -f2)
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


res_file="${folder}/details.csv"

rawFolder="${folder}/raw"
mkdir -p ${rawFolder}

#get machine env
./machine-state.sh > ${rawFolder}/machine-state.txt


ctr=0
for matrix in $matrix_name; do
    raw_file="${rawFolder}/${matrix}.txt"
    tmpFile="${rawFolder}/${matrix}.tmp"
    rcmFlag=""

    KMP_WARNINGS=0 MKL_NUM_THREADS=$thread OMP_NUM_THREADS=${thread} \
        OMP_SCHEDULE=static OMP_PROC_BIND=close OMP_PLACES=cores \
        taskset -c 0-$((thread-1)) ${execFolder}/matStat \
        -m "${matrixFolder}/${matrix}" -c ${thread} -T 1e-4 ${rcmFlag} > ${tmpFile}

    cat ${tmpFile} >> ${raw_file}
    if [[ $printHead == "1" ]]; then
        if [[ ${ctr} == 0 ]]; then
            echo "Matrix,Thread,Nrows,NNZ,symmetric,isDiagZero" > ${res_file}
        fi
    fi

    nrows_w_space=$(cat ${tmpFile} | grep "Nrows =" | cut -d"=" -f2)
    nrows=$(echo ${nrows_w_space})

    nnz_w_space=$(cat ${tmpFile} | grep "NNZ =" | cut -d"=" -f2)
    nnz=$(echo ${nnz_w_space})

    symmetric_w_space=$(cat ${tmpFile} | grep "Symmetric =" | cut -d"=" -f2)
    symmetric=$(echo ${symmetric_w_space})

    isDiagZero_w_space=$(cat ${tmpFile} | grep "is any diag zero =" | cut -d"=" -f2)
    isDiagZero=$(echo ${isDiagZero_w_space})

    echo "${matrix},${thread},${nrows},${nnz},${symmetric},${isDiagZero}" >> ${res_file}
    let ctr=${ctr}+1
    rm -rf ${tmpFile}
done

#store aligned CSV file
sed 's/,/:,/g' ${res_file} | column -t -s: | sed 's/ ,/,/g' > "${alignFolder}/details.csv"

#compress raw folder
cd ${folder}
tar -cvzf raw.tar.gz raw
rm -rf raw
cd -
