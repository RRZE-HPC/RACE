col_types="ABMC MC RACE"

for col in ${col_types}; do
    cp symm_spmv_best_w_likwid.batch.template symm_spmv_best_w_likwid_${col}.batch
    sed -i "s/@COL_TYPE@/${col}/g" symm_spmv_best_w_likwid_${col}.batch
    sbatch symm_spmv_best_w_likwid_${col}.batch
done
