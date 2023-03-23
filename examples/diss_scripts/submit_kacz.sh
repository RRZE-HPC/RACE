col_types="ABMC MC RACE"

for col in ${col_types}; do
    cp kacz.batch.template kacz_${col}.batch
    sed -i "s/@COL_TYPE@/${col}/g" kacz_${col}.batch
    sbatch kacz_${col}.batch
done
