# Variables for running
prc := 4
nx := 48
ny := 48
nz := 96
itrs := 100
dt := 0.002

# Files and flags for the compiler
CC := mpicc
CFLAGS := -O0 -Wall
IFFTW := /m100_work/PROJECTS/spack/spack-0.14/install/linux-rhel7-power9le/gcc-8.4.0/fftw-3.3.8-hwlrarpm6cvjlukhfdowwveb7g7oqwgc/include/
LFFTW := /m100_work/PROJECTS/spack/spack-0.14/install/linux-rhel7-power9le/gcc-8.4.0/fftw-3.3.8-hwlrarpm6cvjlukhfdowwveb7g7oqwgc/lib/ 
INCLUDE := -I./include/ -I$(IFFTW) 
LINK := -L$(LFFTW) -lm -lfftw3
MAIN := $(wildcard *.c)
EXE := $(MAIN:.c=.x)
SRC := $(wildcard src/*.c) $(MAIN)
OBJ := $(SRC:.c=.o)

ifeq ($(flag), debug)
	CFLAGS += -D_DEBUG
endif

# Generating the executables
all: $(EXE)

fftw3_mpi: CFLAGS += -D_FFTW3_MPI
fftw3_mpi: EXE := fftw3_mpi_$(EXE)
fftw3_mpi: LINK += -lfftw3_mpi 
fftw3_mpi: all

# Compiling the object files
%.o: %.c 
	$(CC) -c $< -o $@ $(CFLAGS) $(INCLUDE) 

# Linking the executable
$(EXE): $(OBJ) 
	$(CC) -o $(EXE) $^ $(LINK) $(CFLAGS)
	@rm $(OBJ)

flush:
	@rm -f data/*.dat plots/*.png plots/*.gif

run:
	@mpirun -np $(prc) ./$(EXE) $(itrs) $(nx) $(ny) $(nz) $(dt)

# Running the executable
fftw3_mpi_run: EXE := fftw3_mpi_$(EXE)
fftw3_mpi_run: run

clean: 
	@rm -f $(OBJ) *.x *~
	
plot:
	gnuplot plots/animate.plt
	
format: $(SRC) ./include/utilities.h
	@clang-format -i $^  -verbose || echo "Please install clang-format to run this command"
	
.PHONY: clean plot run flush
