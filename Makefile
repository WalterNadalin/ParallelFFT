# Variables for running
prc := 4
nx := 48
ny := 48
nz := 96
itrs := 101
dt := 0.002

# Files and flags for the compiler
CC := mpicc
CFLAGS := -O0 -Wall
INCLUDE := -I./include/ -I/usr/local/include
LINK := -L/usr/local/lib -lfftw3_mpi -lfftw3 -lm
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
	mpirun -np $(prc) ./$(EXE) $(itrs) $(nx) $(ny) $(nz) $(dt)

# Running the executable
fftw3_mpi_run: EXE := fftw3_mpi_$(EXE)
fftw3_mpi_run: run

clean: 
	@rm -f src/*.o *.x *~
	
plot:
	gnuplot plots/animate.plt
	
format: $(SRC) ./include/utilities.h
	@clang-format -i $^  -verbose || echo "Please install clang-format to run this command"
	
.PHONY: clean plot run flush
