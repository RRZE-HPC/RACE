# README #

RACE-Recursive Algebraic Coloring Engine is a graph coloring library that helps in parallelization of sparse matrix kernels having distance-k (k>=1)
dependencies.RACE uses a recursive level based method for coloring.

### FEATURES ###

* Hardware friendly
* Solve Distance-1 & Distance-2 dependencies kernels
* Supports pre-processing and processing phase
* Easy parallelisation, user needs to just supply serial kernel 
* Self-pinning

### How do I get set up? ###

* git clone the repository
* CC=icc CXX=icpc cmake ..
* Configure the library using ccmake .
* make
* make install
* Library Dependencies : hwloc (will be cloned and installed if not found)


### Contribution guidelines ###



### Who do I talk to? ###

* Christie Louis Alappat <christie.alappat@fau.de>