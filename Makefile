# Variables for running
prc   ?= 4
nx    ?= 32
ny    ?= 32
nz    ?= 96
itrs  ?= 100
dt    ?= 0.001
debug ?= no

# Flags for the compiler
HOSTNAME := $(shell hostname)
CC       := mpicc
CFLAGS   := -O0 -Wall

# Libraries
ifeq ($(HOSTNAME), pop-os) # On my personal computer
	FFTW := usr/local
else # On Marconi100
	FFTW := m100_work/PROJECTS/spack/spack-0.14/install/linux-rhel7-power9le/gcc-8.4.0/fftw-3.3.8-hwlrarpm6cvjlukhfdowwveb7g7oqwgc
endif

INCLUDE := -I./include/ -I/$(FFTW)/include 
LINK := -L/$(FFTW)/lib -lfftw3_mpi -lfftw3 -lm

# Files
MAIN := $(wildcard *.c)
EXE := $(MAIN:.c=.x)
SRC := $(wildcard src/*.c) $(MAIN)
OBJ := $(SRC:.c=.o)

ifeq ($(debug), yes)
	CFLAGS += -D_DEBUG
	EXE := $(EXE:.x=_debug.x)
endif

# Generating the executables
all: $(EXE)

fftw3_mpi: CFLAGS += -D_FFTW3_MPI
fftw3_mpi: LINK   += -lfftw3_mpi 
fftw3_mpi: fftw3_mpi$(EXE)

# Compiling the object files
%.o: %.c 
	$(CC) -c $< -o $@ $(CFLAGS) $(INCLUDE) 

# Linking the executable
%.x: $(OBJ) 
	$(CC) -o $@ $^ $(LINK) $(CFLAGS)
	@rm $(OBJ)

run: $(EXE)
	@mpirun -np $(prc) ./$^ $(itrs) $(nx) $(ny) $(nz) $(dt)

# Running the executable
%run: %
	@mpirun -np $(prc) ./$^$(EXE) $(itrs) $(nx) $(ny) $(nz) $(dt)

clean: 
	@rm -f $(OBJ) *.x 

flush:
	@rm -f data/*concentration*.dat data/diffusivity.dat plots/*.png
	
plot: data/concentration_1.dat
	gnuplot plots/animate.plt
	
data/concentration_1.dat: CFLAGS += -D_DEBUG
data/concentration_1.dat: $(EXE:.x=_debug.x)
	@mpirun -np $(prc) ./$^ $(itrs) $(nx) $(ny) $(nz) $(dt)

format: $(SRC) ./include/utilities.h
	@clang-format -i $^  -verbose || echo "Please install clang-format to run this command"
	
.PHONY: clean plot run flush all fftw3_mpi
.INTERMEDIATE: $(OBJ) data/concentration_1.dat
