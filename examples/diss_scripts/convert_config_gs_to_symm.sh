files=$(ls kacz_config*)
for file in ${files}; do
    sed -i "s@folder=kacz@folder=symm_kacz@g" $file
done
