col_types="ABMC MC RACE"

for col in ${col_types}; do
    cp gs.batch.template gs_${col}.batch
    sed -i "s/@COL_TYPE@/${col}/g" gs_${col}.batch
    sbatch gs_${col}.batch
done
