RACE_INC_DIR=~/race/lib/include/
RACE_LIB_DIR=~/race/lib/lib/RACE

INC_DIR=~/race/examples/include
DEPS = parse.c sparsemat.c mmio.c densemat.c kernels.c
OBJ = parse.o sparsemat.o densemat.o kernels.o mmio.o main.o

CC=icpc

CFLAGS=-std=c++11 -qopenmp -O3 -xHOST -fno-alias -I . -I $(RACE_INC_DIR) -I $(INC_DIR) $(LIKWID_INC) -DLIKWID_PERFMON
LFLAGS=$(shell $(RACE_LIB_DIR)/RACE_LFlags.sh)

%.o: %.c $(DEPS)
		$(CC) -c -o $@ $< $(CFLAGS)

race: $(OBJ)
		$(CC) -o $@ $^ $(CFLAGS) $(LFLAGS)

clean: 
	rm *.o race
