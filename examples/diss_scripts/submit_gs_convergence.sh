col_types="SERIAL ABMC MC RACE"

for col in ${col_types}; do
    cp gs_convergence.batch.template gs_convergence_${col}.batch
    sed -i "s/@COL_TYPE@/${col}/g" gs_convergence_${col}.batch
    sbatch gs_convergence_${col}.batch
done
