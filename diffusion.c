// Created by G.P. Brandino, I. Girotto, R. Gebauer
// Modified by W. Nadalin
// Last revision: April 2023

#include "utilities.h"

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  const double L1 = 10., L2 = 10., L3 = 20.;                            // Dimensions of the system
  const int nx = atoi(argv[2]), ny = atoi(argv[3]), nz = atoi(argv[4]); // Grid size
  const double dt = atof(argv[5]); // Time step for time integration
  const int nstep = atoi(argv[1]); // Number of time steps
  const double rad_diff = 0.7;     // Radius of diffusion channel
  const double rad_conc = 0.6;     // Radius of starting concentration
  double start, end; // For time measures
  double fxconc, fyconc, fzconc, fxdiff, fydiff, fzdiff, fac, ss, ss_all, xx, xy, xz;
  double *diffusivity, *concentration, *dconc, *auxx, *auxy;
  int size, i, ix, iy, iz, ipol, istep, index;
  fftw_dist_handler fft_h;

  MPI_Comm_size(MPI_COMM_WORLD, &size); // Initialization of the MPI environment
  init_fftw(&fft_h, nx, ny, nz, MPI_COMM_WORLD); // Initialize the fftw system and local dimension

  const int nx_local = fft_h.local_nx, nx_local_offset = fft_h.local_nx_offset,
            global_size_grid = fft_h.global_size_grid, local_size_grid = fft_h.local_size_grid;

  // Allocate distributed memory arrays
  diffusivity = (double *)malloc(local_size_grid * sizeof(double));
  concentration = (double *)malloc(local_size_grid * sizeof(double));
  dconc = (double *)malloc(local_size_grid * sizeof(double));
  auxx = (double *)malloc(local_size_grid * sizeof(double));
  auxy = (double *)malloc(local_size_grid * sizeof(double));

  ss = 0.0; // To integrate (and normalize) the concentration

  for (iz = 0; iz < nz; ++iz) { // Initialization of diffusion coefficent and concentration
    xz = L3 * ((double)iz) / nz;
    fzdiff = exp(-pow((xz - 0.5 * L3) / rad_diff, 2));
    fzconc = exp(-pow((xz - 0.5 * L3) / rad_conc, 2));

    for (iy = 0; iy < ny; ++iy) {
      xy = L2 * ((double)iy) / ny;
      fydiff = exp(-pow((xy - 0.5 * L2) / rad_diff, 2));
      fyconc = exp(-pow((xy - 0.5 * L2) / rad_conc, 2));

      for (ix = 0; ix < nx_local; ++ix) {
        xx = L1 * ((double)(ix + nx_local_offset)) / nx;
        fxdiff = exp(-pow((xx - 0.5 * L1) / rad_diff, 2));
        fxconc = exp(-pow((xx - 0.5 * L1) / rad_conc, 2));

        index = index_f(ix, iy, iz, ny, nz);
        diffusivity[index] = MAX(fxdiff * fydiff, fydiff * fzdiff);
        concentration[index] = fxconc * fyconc * fzconc;
        ss += concentration[index];
      }
    }
  }

  // Normalize the concentration
  fac = L1 * L2 * L3 / global_size_grid;
  MPI_Allreduce(&ss, &ss_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  ss_all = 1.0 / (ss_all * fac);

  for (ix = 0; ix < local_size_grid; ++ix)
    concentration[ix] *= ss_all;

#ifdef _DEBUG //  Printing initial concentration and coefficent
	int how_many = 10;
	
  plot_data_2d("data/diffusivity", nx, ny, nz, nx_local, nx_local_offset, 2, diffusivity);
  plot_data_2d("data/initial_concentration", nx, ny, nz, nx_local, nx_local_offset, 2,
               concentration);
#endif	

  start = seconds();

  for (i = 0; i < local_size_grid; ++i)
    dconc[ix] = 0.0;

  for (istep = 1; istep <= nstep; ++istep) { // Start the dynamics
    for (ipol = 1; ipol <= 3; ++ipol) { // Compute derivative
      derivative(&fft_h, nx, ny, nz, L1, L2, L3, ipol, concentration, auxx);

      for (i = 0; i < local_size_grid; ++i)
        auxx[i] *= diffusivity[i];

      derivative(&fft_h, nx, ny, nz, L1, L2, L3, ipol, auxx, auxy);

      for (i = 0; i < local_size_grid; ++i)
        dconc[i] += auxy[i]; // Summing up contributions from the three spatial directions
    }

    for (i = 0; i < local_size_grid; ++i) { // Update concentration
      concentration[i] += dt * dconc[i];
      dconc[i] = 0.0;
    }

    end = seconds();

#ifdef _DEBUG
    if (istep % how_many == 1) // Check and save data
      print_info(concentration, nx, ny, nz, nx_local, nx_local_offset, L1, L2, L3, istep, 
      					 how_many, end - start);
#endif
  }

#ifndef _DEBUG	
  double time, max_time;
  int rank;
  
  time = end - start;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Allreduce(&time, &max_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);	

  if(!rank) {
#ifdef _FFTW3_MPI
    const char *mode = "fftw_mpi";
#else
    const char *mode = "homemade";
#endif  
		
    const char *times = "data/times.dat"; // Where to write time
    FILE *file = fopen(times, "a");
		
    fprintf(file, "%s\t%d\t%d\t%d\t%d\t%d\t%lf\t%lf\n", mode, size, nx, ny, nz, nstep, dt, max_time);

    fclose(file);
  }
#endif

  close_fftw(&fft_h);
  free(diffusivity);
  free(concentration);
  free(dconc);
  free(auxx);
  free(auxy);

  MPI_Finalize();

  return 0;
}
