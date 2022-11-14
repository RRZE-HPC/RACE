# README #

RACE-Recursive Algebraic Coloring Engine is a graph library serving two main functionalities:

* coloring
* temporal blocking

The coloring part of the library finds its use in parallelization of sparse matrix kernels having distance-k (k>=1) dependencies. For example kernels like Gauss-Seidel, SpMTV, SymmSpMV and Kaczmarz can be parallelized using RACE.
The temporal blocking part allows for cache blocking of multiple calls of a sparse kernel. For example the matrix power kernel (MPK) where multiple SpMVs are used.
RACE works on the graph of the matrix and is an algebraic approach.
Several optimizations have been implemented in RACE to boost the runtime performance of the kernels on modern CPUs.

### FEATURES ###

* Hardware friendly
* Solve Distance-1, Distance-2 dependencies in kernels and allows to cache block MPK-style kernels
* Supports pre-processing and processing phase
* Easy parallelisation and blocking, user needs to just supply serial kernel
* Support for CRS data format

### How to build RACE? ###

* `git clone git@bitbucket.org:essex/race.git`
* `cd race && mkdir build && cd build`
* `CC=$(C\_COMPILER) CXX=$(CXX\_COMPILER) cmake .. -DCMAKE_INSTALL_PREFIX=(where to install)`
* Configure the library using `ccmake .` (if needed)
* `make`
* `make install`
* Library Dependencies : [hwloc](https://www.open-mpi.org/projects/hwloc/), [Intel SpMP](https://github.com/IntelLabs/SpMP), [GAP](https://github.com/sbeamer/gapbs). The libraries will be cloned and installed if needed.
* When using RACE library, use CMAKE find\_package to get the proper linking flags for the library.

### Want to try RACE? ###
RACE provides examples to illustrate the usage and easiness of using the RACE library. 
Please see the [example folder](https://bitbucket.org/essex/race/src/master/examples/).

### Working of RACE coloring ###
![Screenshot](animations/domain_anim.gif)
![Scrrentshot](animations/zone_tree_anim.gif)


### Citing RACE ###

If you are using RACE for coloring please use the following reference:

* C. L. Alappat, G. Hager, O. Schenk, J. Thies, A. Basermann, A. R. Bishop, H. Fehske, and G. Wellein:
  A Recursive Algebraic Coloring Technique for Hardware-Efficient Symmetric Sparse Matrix-Vector Multiplication, published in ACM TOPC journal.  [doi:10.1145/3399732](https://doi.org/10.1145/3399732)

If you are using RACE for matrix power or temporal blocking of sparse matrices please use the following reference:

* C. L. Alappat, G. Hager, O. Schenk, and G. Wellein:
  Level-based Blocking for Sparse Matrices: Sparse Matrix-Power-Vector Multiplication. Accepted for publication in IEEE TPDS. Preprint: [arXiv:2205.01598](https://arxiv.org/abs/2205.01598)


### Further Reading ###

* A poster on RACE was presented at SC'18. The poster and abstract can be found [here](https://sc18.supercomputing.org/proceedings/src_poster/src_poster_pages/spost109.html).
* A short paper on RACE showing how it can be used to parallelise an eigen value solver can be found [here](https://src.acm.org/binaries/content/assets/src/2019/christie-louis-alappat.pdf).
  The paper won the second place in [ACM SRC](https://src.acm.org) 2019.
