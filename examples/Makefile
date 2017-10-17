RACE_INC_DIR=~/race/lib/include/
RACE_LIB_DIR=~/race/lib/lib/RACE
RACE_LIB=

HWLOC_INC_DIR=~/race/build_RACE/hwloc-inst/include
HWLOC_LIB_DIR=~/race/build_RACE/hwloc-inst/lib
HWLOC_LIB=


INC_DIR=~/race/examples/include
DEPS = parse.c sparsemat.c mmio.c densemat.c kernels.c
OBJ = parse.o sparsemat.o densemat.o kernels.o mmio.o main.o

CC=icpc

CFLAGS=-std=c++11 -qopenmp -I . -I $(RACE_INC_DIR) -I $(INC_DIR) $(LIKWID_INC) -qopenmp -DLIKWID_PERFMON $(RACE_LIB_DIR)/libRACE.a $(LIKWID_LIB)

ifeq ($(STATIC),true)
	CFLAGS=-std=c++11 -qopenmp -I . -I $(RACE_INC_DIR) -I $(INC_DIR) $(LIKWID_INC) -qopenmp -DLIKWID_PERFMON $(RACE_LIB_DIR)/libRACE.a $(LIKWID_LIB)
else
	CFLAGS=-std=c++11 -qopenmp -I . -I $(RACE_INC_DIR) -I $(INC_DIR) $(LIKWID_INC) -qopenmp -DLIKWID_PERFMON $(RACE_LIB_DIR)/libRACE.a $(LIKWID_LIB)
CFLAGS=-std=c++11 -qopenmp -I . -I $(RACE_INC_DIR) -I $(INC_DIR) $(LIKWID_INC) -qopenmp -DLIKWID_PERFMON -L $(RACE_LIB_DIR) -L $(HWLOC_LIB_DIR) -Wl,-rpath $(HWLOC_LIB_DIR) -Wl,-rpath $(RACE_LIB_DIR) $(LIKWID_LIB)
	RACE_LIB=-l RACE -l hwloc
endif


%.o: %.c $(DEPS)
		$(CC) -c -o $@ $< $(CFLAGS)

race: $(OBJ)
		$(CC) -o $@ $^ $(CFLAGS) $(RACE_LIB)  -l likwid

clean: 
	rm *.o race
