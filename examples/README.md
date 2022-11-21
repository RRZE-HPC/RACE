### How to compile RACE examples? ###
RACE provides examples to illustrate the usage and easiness of using the RACE library. To try it out:

* `cd race/example`
* `mkdir build && cd build`
* `CC=$(C\_COMPILER) CXX=$(CXX\_COMPILER) cmake .. -DRACE\_DIR=$(RACE\_LIB)`
* optionally set `RACE_USE_SPMP` for pre-permuting matrix with RCM permutation.
  [Intel SpMP](https://github.com/IntelLabs/SpMP) library is used for RCM
  permutation.
* `make`

Note that as some examples depend on sparse BLAS library, linking with BLAS is necessary. The CMake tries to automatically detect BLAS library through environment flags, for example MKLROOT of Intel MKL. However, if RACE couldn't find a BLAS library please link it manually by specifying the include directory using `-DCBLAS_INCLUDE_DIR=<path to BLAS include directory>`.

### Running coloring examples ###

* To run benchmarks SpMTV, SymmSpMV, KACZ, GS with 10 threads on first 10 cores, and &epsilon;<sub>0</sub>= &epsilon;<sub>1</sub>=0.8
  use the following command:
  ```bash
  RACE_EFFICIENCY=80,80 taskset -c 0-9 ./coloring  -m <matrix market file> -c 10 -t 1
  ```
* Note pinning is carried out by RACE in this case
* To get other options use: ./coloring -h
* Relevant environment variables are:

  | Name | Values | Default | Recomendation | Purpose | Example
  | :---   | :--- | :--- | :--- | :--- | :--- |
  | `RACE_EFFICIENCY` | [40, 100] | 40 | If you are unsure, set 80 for first few stages  | Controls recursion. It sets the load imbalance tolerance &epsilon;<sub>s</sub> at each recursive stage s | `RACE_EFFICIENCY=50,80` implies &epsilon;<sub>0</sub>=0.5 and &epsilon;<sub>1</sub>=0.8.
  | `RACE_THREADS` | [1, c], where c is the total number of cores. | auto | If you are unsure, don't specify this variable. | Assign n number of threads at a particular recusion stage. This overrides the internal heuristic to determine number of threads at each stage.| `RACE_THREADS=2,3` implies 2 threads in stage 0 and 3 threads in stage 1 of recursion.
  
### Running temporal blocking examples ###

* To run MPK kernel with power 3, i.e., A<sup>3</sup>x computation, with 10 threads on first 10 cores on an architecture with last level cache size of 40 MB, 
  use the following command:
  ```bash
  RACE_CACHE_VIOLATION_CUTOFF=1,1 OMP_NUM_THREADS=10 OMP_PLACES=cores OMP_PROC_BIND=close  taskset -c 0-9 ./mtxPower -m <matrix market file> -s 40 -i 3
  ```
* Note in this case pinning is left to user, therefore we use here OMP environment variables
* To get other options use: ./mtxPower -h
* Relevant environment variables are:

  | Name | Values | Default | Recomendation | Purpose | Example
  | :---   | :--- | :--- | :--- | :--- | :--- |
  | `RACE_CACHE_VIOLATION_CUTOFF` | [1,inf] | see next | If you are unsure, set 1 for first few stages  | Controls recursion. It sets the safety factor of cache at different recursion stage. Higher the value less the probablity of going for further recursion. | `RACE_CACHE_VIOLATION_CUTOFF=1,2` implies safety factor of 1 at stage 0 and 2 at stage 1.
  | `RACE_CACHE_VIOLATION_CUTOFF_DEFAULT` | [1, inf] | 50 |If you are unsure, don't specify this variable. | Sets the default safety factor value of all stages. | `RACE_CACHE_VIOLATION_CUTOFF_DEFAULT=2` implies safety factor is 2 by default on all stages. 
  | `RACE_MAX_RECURSION_STAGES` | [1, inf] | auto |If you are unsure, don't specify this variable. | Sets the maximum number of recursion stages. | `RACE_MAX_RECURSION_STAGES=1` implies only 1 recursion stage is allowed.
