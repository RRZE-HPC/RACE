matrixFolder="/home/vault/unrz/unrz002h/matrix/RACE_test"

#Checks performance of some matrices, once some new changes are made
./check-state.sh config.txt

#Step 1: Check if there are some new commits
lc=$(cat last_commit.txt | head -n 1)
cc=$(git log --pretty=format:"%h" | head -n 1)

if [[ "$lc" == "$cc" ]]; then
    echo "There is no new commit, so no new measurement"
    exit 2
else
    echo $cc > last_commit.txt
    git log --pretty=oneline --abbrev-commit | head -n 1 >> commit_description.txt
    git log --pretty=format:"%h %ad" --abbrev-commit | head -n 1 >> commit_date.txt
fi

echo $cc
folder="${cc}"
mkdir -p $folder
mkdir -p $folder/raw

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

#Only full socket run
threads=$(cat $folder/machine-state.txt | grep "Cores per socket" | cut -d":" -f2)
#powers to test
powers="1 2 3 4 6 8 10"

race_file="${folder}/race.txt"
printf "%20s, %10s, %10s, %10s, %10s, %10s\n"\
    "Matrix" "SpMV" "SpMTV" "GS" "KACZ" "SymmSpMV" > ${race_file}

mtxPower_file="${folder}/mtxPower.txt"
printf "%20s" "Matrix" > ${mtxPower_file}
for power in ${powers}; do
    printf ", %10s" "p=$power" >> ${mtxPower_file}
done
printf "\n" >> ${mtxPower_file}

function readResult
{
    cat $1 | grep "$2" | head -n 1 | cut -d "=" -f 2 | cut -d "G" -f1
}

for matrix in $matrix_name; do
    raw_file="${folder}/raw/${matrix}.txt"
    RACE_EFFICIENCY=80,80 OMP_NUM_THREADS=$threads taskset -c 0-$((threads-1)) ./race -m ${matrixFolder}/${matrix} \
        -c $threads -t 1 -v > tmp.txt
    cat tmp.txt > ${raw_file}

    spmv=$(readResult "tmp.txt" "SPMV")
    spmtv=$(readResult "tmp.txt" "SPMTV")
    gs=$(readResult "tmp.txt" "GS")
    kacz=$(readResult "tmp.txt" "KACZ")
    symm_spmv=$(readResult "tmp.txt" "SYMM_SPMV")
    printf "%20s, %10.4f, %10.4f, %10.4f, %10.4f, %10.4f\n"\
        "${matrix}" $spmv $spmtv $gs $kacz $symm_spmv >> ${race_file}

    printf "%20s" "${matrix}" >> ${mtxPower_file}
    for power in $powers; do
        OMP_NUM_THREADS=$threads taskset -c 0-$((threads-1)) ./mtxPower \
            -m ${matrixFolder}/${matrix} -c $threads -t 1 \
            -i $power -n 1 -v > tmp.txt
        cat tmp.txt >> ${raw_file}
        perf=$(readResult "tmp.txt" "RACE power perf.")
        printf ", %10.4f" ${perf} >> ${mtxPower_file}
    done
    printf "\n" >> ${mtxPower_file}
done

rm -rf tmp.txt

#now generate plot

#Modify config file
max_commits=4
legends=$(cat "commit_date.txt" | head -n ${max_commits} | cut -d" " -f1)
rest_cols=""
latest_commit=$(cat "commit_date.txt" | head -n 1 | cut -d" " -f1)
spmv_file="../${latest_commit}/race.txt"
commit_files=""
commit_files_power=""
for legend in $legends; do
    rest_cols="${rest_cols}6 " #SymmSpMV is 6-th column
    commit_files="${commit_files}../${legend}/race.txt "
    commit_files_power="${commit_files_power}../${legend}/mtxPower.txt "
done

#plot symmspmv
cd plot_generator
cp config_symmspmv.txt curr_config.txt
./substitute.sh "spmv_col" "2" curr_config.txt
./substitute.sh "rest_cols" "${rest_cols}" curr_config.txt
./substitute.sh "spmv_file" "${spmv_file}" curr_config.txt
./substitute.sh "commit_files" "${commit_files}" curr_config.txt
./substitute.sh "rest_legends" "${legends}" curr_config.txt
cat curr_config.txt
./generateRACE_plot.sh curr_config.txt
rm curr_config.txt

#plot power
rm -rf powerPlots
mkdir -p powerPlots
curr_ctr=2
#remember data file will be transposed by the plot script,
#so column numbers should correspond to that
for matrix in ${matrix_name}; do
    cp config_mtxpower.txt curr_config.txt
    curr_cols=""
    for file in ${commit_files_power}; do
        curr_cols="${curr_cols}${curr_ctr} "
    done
    ./substitute.sh "cols" "$curr_cols" curr_config.txt
    ./substitute.sh "commit_files" "${commit_files_power}" curr_config.txt
    ./substitute.sh "legends" "${legends}" curr_config.txt
    escaped_matrix_name=$(sed -e "s/_/-/g" <<< $matrix)
    ./substitute.sh "matrix" "${escaped_matrix_name}" curr_config.txt
    cat curr_config.txt
    ./generatePower_plot.sh curr_config.txt
    rm curr_config.txt
    let curr_ctr=${curr_ctr}+1
done

./generatePower_plot_all.sh
cd -
