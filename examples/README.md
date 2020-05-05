### How to run RACE examples? ###
RACE provides examples to illustrate the usage and easiness of using the RACE library. To try it out:

* `cd race/example`
* `mkdir build && cd build`
* `CC=$(C\_COMPILER) CXX=$(CXX\_COMPILER) cmake .. -DRACE\_DIR=$(RACE\_LIB)`
* optionally set `RACE_USE_SPMP` for pre-permuting matrix with RCM permutation.
  [Intel SpMP](https://github.com/IntelLabs/SpMP) library is used for RCM
  permutation.
* `make`
* To run benchmarks SpMTV, SymmSpMV, KACZ, GS with 10 threads on first 10 cores, and \epsilon_0=\epsilon_1=0.8 use the following command:

```
RACE_EFFICIENCY=80,80 taskset -c 0-9 ./race  -m <matrix market file> -c 10 -t 1
```
* To get other options use: ./race -h
