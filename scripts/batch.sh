#!/bin/bash
#SBATCH -A tra23_units
#SBATCH -p m100_usr_prod
#SBATCH --time 02:00:00       # format: HH:MM:SS
#SBATCH -N 8                  # nodes
#SBATCH --gres=gpu:4          # gpus per node out of 4
#SBATCH --mem=246000          # memory per node out of 246000MB
#SBATCH --ntasks-per-node=32  # 8 tasks out of 128
#SBATCH --ntasks-per-core=1
#SBATCH --job-name=wanda_parallel_jacobi
#SBATCH --mail-type=ALL
#SBATCH --mail-user=walter.nadalin@studenti.units.it

module purge
module load spectrum_mpi

make clean
make 
make fftw3_mpi

rm data/times.dat
echo -e "mode\t\t\tprc\tnx\tny\tnz\titr\tdt\t\ttime" >> data/times.dat

nx=48
ny=48
nz=96
nstep=101
dt=0.002
prc=32

for value in {1..4}
do
        mpirun -np $value ./diffusion.x $nstep $nx $ny $nz $dt #-npernode 32
        ((prc*=2))
done

for value in {1..4}
do
        mpirun -np $value ./fftw3_mpi_diffusion.x $nstep $nx $ny $nz $dt
        ((prc*=2))
done

make clean

