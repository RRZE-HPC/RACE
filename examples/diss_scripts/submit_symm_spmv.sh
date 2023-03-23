col_types="ABMC MC RACE"

for col in ${col_types}; do
    cp symm_spmv.batch.template symm_spmv_${col}.batch
    sed -i "s/@COL_TYPE@/${col}/g" symm_spmv_${col}.batch
    sbatch symm_spmv_${col}.batch
done
