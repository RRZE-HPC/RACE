GAP_folder=$1
fileToPatch=${GAP_folder}/src/platform_atomics.h
cp ${fileToPatch} ${GAP_folder}/src/platform_atomics.h.unpatched
toPatchLine=$(cat ${fileToPatch} | grep -n "bool compare_and_swap(float &x," | cut -d":" -f1 | head -n1)
sed -i "$((toPatchLine-1))s@template<>@template<> inline@g" ${fileToPatch}
toPatchLine=$(cat ${fileToPatch} | grep -n "bool compare_and_swap(double &x," | cut -d":" -f1 | head -n1)
sed -i "$((toPatchLine-1))s@template<>@template<> inline@g" ${fileToPatch}
