col_types="SERIAL ABMC MC RACE"

for col in ${col_types}; do
    cp kacz_convergence.batch.template kacz_convergence_${col}.batch
    sed -i "s/@COL_TYPE@/${col}/g" kacz_convergence_${col}.batch
    sbatch kacz_convergence_${col}.batch
done
