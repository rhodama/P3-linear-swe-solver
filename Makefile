CPP=CC
CFLAGS= -O3 -lm -march=native -funroll-loops -mavx2
OPTFLAGS=-O3  
MPIFLAGS=-DMPI 
DEBUGFLAGS=-g -pg

NVCC=nvcc
NVCCFLAGS= -Xptxas -dlcm=ca -O3 -DCUDA -gencode arch=compute_80,code=sm_80 
PYTHON=python3

all: mpi gpu basic_serial serial

mpi: build/mpi
gpu: build/gpu
serial: build/serial
basic_serial: build/basic_serial

build/mpi: common/main.cpp common/scenarios.cpp mpi/mpi.cpp
	$(CPP) $^ -o $@ $(MPIFLAGS) $(CFLAGS) $(OPTFLAGS)

build/gpu: common/main.cpp common/scenarios.cpp gpu/gpu.cu
	$(NVCC) $^ -o $@ $(NVCCFLAGS)

build/serial: common/main.cpp common/scenarios.cpp serial/serial.cpp
	$(CPP) $^ -o $@ $(CFLAGS) $(COPTFLAGS)

build/basic_serial: common/main.cpp common/scenarios.cpp serial/basic_serial.cpp
	$(CPP) $^ -o $@ $(CFLAGS) $(COPTFLAGS)

.PHONY: clean

clean:
	rm -f build/*.out
	rm -f build/*.o
	rm -f build/*.gif