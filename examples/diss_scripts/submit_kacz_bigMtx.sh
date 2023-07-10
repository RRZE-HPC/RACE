col_types="ABMC MC RACE"

for col in ${col_types}; do
    cp kacz_bigMtx.batch.template kacz_${col}_bigMtx.batch
    sed -i "s/@COL_TYPE@/${col}/g" kacz_${col}_bigMtx.batch
    sbatch kacz_${col}_bigMtx.batch
done
