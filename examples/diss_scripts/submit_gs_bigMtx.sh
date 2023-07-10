col_types="ABMC MC RACE"

for col in ${col_types}; do
    cp gs_bigMtx.batch.template gs_${col}_bigMtx.batch
    sed -i "s/@COL_TYPE@/${col}/g" gs_${col}_bigMtx.batch
    sbatch gs_${col}_bigMtx.batch
done
