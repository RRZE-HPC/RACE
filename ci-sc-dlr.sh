#!/usr/bin/env bash
set -e

## default options and declarations
# kernel lib
PRGENV="gcc-4.9.2-openmpi-1.10.1" # intel-13.0.1-mpich gcc-4.8.2-openmpi
BUILD_TYPE=Release
INSTALL_PREFIX=../../

# list of modules to load
MODULES_BASIC="cmake ccache"

## parse command line arguments
usage() { echo "Usage: $0 [-e <PrgEnv/module-string>] [-b <Release|Debug|...>]" 1>&2; 
exit 1; }

while getopts "e:b:p:h" o; do
    case "${o}" in
        e)
            PRGENV=${OPTARG}
            ;;
        b)
            BUILD_TYPE=${OPTARG}
            ;;
        p)
            INSTALL_PREFIX=${OPTARG}
            ;;
        h)
            usage
            ;;
        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))

echo "Options: PRGENV=${PRGENV}, BUILD_TYPE=${BUILD_TYPE}"

## prepare system for compilation
# configure modulesystem
export MODULEPATH=/tools/modulesystem/modulefiles
module() { eval `/usr/bin/modulecmd bash $*`; }

# load modules
module load "PrgEnv/$PRGENV"
# set compiler names
if [[ "$PRGENV" =~ gcc* ]]; then
    export FC="gfortran" CC="gcc" CXX="g++"
elif [[ "$PRGENV" =~ intel* ]]; then
  export FC=ifort CC=icc CXX=icpc
else
  set -- $(mpicc -show)
  export CC=$1
  set -- $(mpicxx -show)
  export CXX=$1
  set -- $(mpif90 -show)
  export FC=$1
fi

echo "compilers: CC=$CC, CXX=$CXX, FC=$FC"

for m in $MODULES_BASIC; do module load $m; done

module list


# "gcc -fsanitize=address" requires this
ulimit -v unlimited


build_dir=build_${PRGENV}_${BUILD_TYPE}
build_examples_dir=build_examples_${PRGENV}_${BUILD_TYPE}
install_dir=$INSTALL_PREFIX/install-${PRGENV}-${BUILD_TYPE}

error=0

function update_error { 
if [[ "${error}" = "0" ]]; then
  error=$1
fi
}


# build and install
mkdir $build_dir                          || update_error ${LINENO}
cd $build_dir                             || update_error ${LINENO}
cmake -DCMAKE_INSTALL_PREFIX=$install_dir \
-DCMAKE_BUILD_TYPE=${BUILD_TYPE} -DBUILD_SHARED_LIBS=ON ..              || update_error ${LINENO}

make -j 24 || make || update_error ${LINENO}
make install       || update_error ${LINENO}
cd ..              || update_error ${LINENO}

# build examples
mkdir ${build_examples_dir} || update_error ${LINENO}
cd ${build_examples_dir} || update_error ${LINENO}
export CMAKE_PREFIX_PATH=$install_dir:${CMAKE_PREFIX_PATH}
cmake ../examples || update_error ${LINENO}
make -j || make    || update_error ${LINENO}
gunzip -k -c ../examples/matrices/spinSZ12.mm.gz > spinSZ12.mm || update_error ${LINENO}

# a simple test, just run with some example matrix
./race -v -m spinSZ12.mm -c 12 -i 25 &> test.log || update_error ${LINENO}
grep -i "validated coloring" test.log || update_error ${LINENO}


#return error code
exit $error
