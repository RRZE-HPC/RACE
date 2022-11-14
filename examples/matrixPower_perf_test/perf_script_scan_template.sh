############# Settings ###################

num=$1

matrixFolder="/home/vault/unrz/unrz002h/matrix/MatrixPower_bigMtx"
#folder="mtxPower_results_Flan_rec_${num}"
folder="mtxPower_results_bigMtx"
tmpFile="tmp_${num}.txt"
powers="1 2 3 4 6"
#powers="4"
threads="20"
#cacheSizes="1 20 25 35 45 55"
#cacheSizes="1 25 35 45 55 65 75 85 95 105"
cacheSizes="35 45"
#safetyFactors="1000 1,1000 1,1,1000 1,1,1,1,1000 1,1,1,1,1,1,1000 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,10000"
#safetyFactors="100000"
stages="1 80"
RCMs="0 1"
trad="off"
node=1

############ Don't change anything below ############################

#./check-state.sh config.txt

mkdir -p $folder
mkdir -p $folder/raw_tradPower
mkdir -p $folder/raw_mtxPower

#get machine env
./machine-state.sh > $folder/machine-state_old.txt
MachineState/machinestate.py -o $folder/machine-state.json

#get matrix names
cd $matrixFolder
matrix_name=

while read -r i
do
    matrix_name=$matrix_name" "$i
done < <(find *)

echo $matrix_name

cd -


function readResult
{
    file=$1
    str=$2
    place=$3
    cutStr=$4
    power=$5
    cs=$6

    if [[ "$cs" == "-1" ]]; then
        cat $file | grep "$str" | grep "power = $power" | head -n 1 | cut -d "=" -f $place| cut -d "${cutStr}" -f1
    else
        cat $file | grep "$str" | grep "cache size = $cs" | grep "power = $power" | head -n 1 | cut -d "=" -f $place | cut -d "${cutStr}" -f1
    fi
}

if [[ ${trad} == "on" ]]; then
    tradPower_file="${folder}/tradPower.csv"
    printf "%20s, %10s, %5s, %5s, %5s, %10s\n" "# Matrix" "Perf" "Thread" "Power" "RCM" "Total_Pre-time" > ${tradPower_file}

    for matrix in $matrix_name; do
        raw_file="${folder}/raw_tradPower/${matrix}.txt"
        for thread in $threads; do
            for RCM in $RCMs; do
                rcmFlag=""
                if [[ $RCM == "1" ]]; then
                    rcmFlag="-R"
                fi

                powers_w_comma=$(echo ${powers} | sed -e "s@ @,@g")
                OMP_NUM_THREADS=$thread likwid-pin -c 0-$((thread-1)) ../tradPower_scan \
                    -m ${matrixFolder}/${matrix} -c $thread -t 1 \
                    -i $powers_w_comma -v $rcmFlag > ${tmpFile}

                cat ${tmpFile} >> ${raw_file}

                for power in ${powers}; do
                    perf=$(readResult "${tmpFile}" "Trad. perf summary:" 3 "G" ${power} -1)
                    preTime=$(readResult "${tmpFile}" "Total pre-processing time" 3 "s" ${power} -1)
                    printf "%20s, %10.4f, %5d, %5d, %5d, %10.6f\n" ${matrix} ${perf} ${thread} ${power} ${RCM} ${preTime} >> ${tradPower_file}
                done
            done
        done
    done
fi


mtxPower_file="${folder}/mtxPower.csv"
printf "%20s, %10s, %5s, %5s, %5s, %5s, %5s, %5s, %10s, %10s\n" "# Matrix" "Perf" "Thread" "Power" "CS" "N" "Stage" "RCM" "RACE_Pre-time" "Total_Pre-time" > ${mtxPower_file}


for matrix in $matrix_name; do
    raw_file="${folder}/raw_mtxPower/${matrix}.txt"
    for thread in $threads; do
        for stage in $stages; do
            safetyFactor=""

            for((idx=0; idx<${stage}; ++idx)); do
                safetyFactor="${safetyFactor}1,"
            done
            safetyFactor="${safetyFactor}100000"

            for RCM in $RCMs; do
                rcmFlag=""
                if [[ $RCM == "1" ]]; then
                    rcmFlag="-R"
                fi

                powers_w_comma=$(echo ${powers} | sed -e "s@ @,@g")
                cacheSize_w_comma=$(echo ${cacheSizes} | sed -e "s@ @,@g")
                RACE_CACHE_VIOLATION_CUTOFF=${safetyFactor} OMP_NUM_THREADS=$thread likwid-pin -c 0-$((thread-1)) ../mtxPower_scan \
                    -m ${matrixFolder}/${matrix} -c $thread -t 1 \
                    -i $powers_w_comma -s $cacheSize_w_comma $rcmFlag > ${tmpFile} #-v

                cat ${tmpFile} >> ${raw_file}

                for power in ${powers}; do
                    for cacheSize in ${cacheSizes}; do
                        perf=$(readResult "${tmpFile}" "RACE perf summary:" 4 "G" ${power} ${cacheSize})
                        #safetyFactor_wo_comma=$(echo ${safetyFactor} | sed -e "s@,@;@g")
                        racePreTime=$(readResult "${tmpFile}" "RACE pre-processing time" 4 "s" ${power} ${cacheSize})
                        preTime=$(readResult "${tmpFile}" "Total pre-processing time" 4 "s" ${power} ${cacheSize})
                        printf "%20s, %10.4f, %5d, %5d, %5.2f, %5d, %5d, %5d, %10.6f, %10.6f\n" ${matrix} ${perf} ${thread} ${power} ${cacheSize} ${node} ${stage} ${RCM} ${racePreTime} ${preTime} >> ${mtxPower_file}
                    done
                done

            done
        done
    done
done

rm -rf ${tmpFile}

cd ${folder}
tar -cvzf raw_mtxPower.tar.gz raw_mtxPower
tar -cvzf raw_tradPower.tar.gz raw_tradPower
rm -rf raw_mtxPower
rm -rf raw_tradPower
cd -
