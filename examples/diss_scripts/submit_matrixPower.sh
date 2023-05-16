nums="1 2 3"

for num in ${nums}; do
    cp matrixPower.batch.template matrixPower_${num}.batch
    sed -i "s/@NUM@/${num}/g" matrixPower_${num}.batch
    sbatch matrixPower_${num}.batch
done
