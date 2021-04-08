urlFile="matrixURL.txt"
  
matrixNames=$(cat ${urlFile} | cut -d"," -f1)
urls=$(cat ${urlFile} | cut -d"," -f2)

matCtr=1
for mat in ${matrixNames}; do
        matURL=$(echo ${urls} | cut -d" " -f${matCtr})
        wget ${matURL}
        tar -xvf "${mat}.tar.gz"
        mv ${mat}/${mat}.mtx .
        rm -rf "${mat}" "${mat}.tar.gz"
        let matCtr=${matCtr}+1
done
