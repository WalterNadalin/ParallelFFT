// Created by G.P. Brandino, I. Girotto, R. Gebauer
// Modified by W. Nadalin
// Last revision: April 2023

#include "utilities.h"

void print_info(double *concentration, int nx, int ny, int nz, int nx_local, int nx_local_offset,
                double L1, double L2, double L3, int istep, double start, double end) {
  double ss = 0., r2mean = 0., ss_all, r2mean_all;
  double fac = L1 * L2 * L3 / (nx * ny * nz), rr, xx, xy, xz;
  int iz, ix, iy, rank, index;

  for (iz = 0; iz < nz; ++iz) { // Check the normalization of concentration
    xz = L3 * ((double)iz) / nz - 0.5 * L3;
    for (iy = 0; iy < ny; ++iy) {
      xy = L2 * ((double)iy) / ny - 0.5 * L2;
      for (ix = 0; ix < nx_local; ++ix) {
        xx = L1 * ((double)(ix + nx_local_offset)) / nx - 0.5 * L1;
        rr = pow(xx, 2) + pow(xy, 2) + pow(xz, 2);
        index = index_f(ix, iy, iz, ny, nz);
        ss += concentration[index];
        r2mean += concentration[index] * rr;
      }
    }
  }

  MPI_Allreduce(&ss, &ss_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&r2mean, &r2mean_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  ss_all = fac * ss_all;
  r2mean_all = fac * r2mean_all;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (!rank)
    printf("%d\t%17.15f %17.15f Elapsed time per iteration %f \n", istep, r2mean_all, ss_all,
           (end - start) / istep);

  // char title[80]; // Title of the plot
  // sprintf(title, "data/concentration_%d", 1 + (istep - 1) / 30);
  // plot_data_2d(title, nx, ny, nz, nx_local, nx_local_offset, 2, concentration);
}

void plot_data_2d(char *name, int nx, int ny, int nz, int nx_local, int nx_local_offset, int dir,
                  double *data) {
  int ix, iy, iz, i;
  FILE *fp;
  char buf[256];
  int index;

  int mype, npes, owner;
  int *sizes, *displ;
  double *buffer, *buffer1d, *local_buffer;

  MPI_Comm_rank(MPI_COMM_WORLD, &mype);
  MPI_Comm_size(MPI_COMM_WORLD, &npes);

  snprintf(buf, sizeof(buf), "%s.dat", name);

  owner = npes + 1;

  // Finding the owner of the middle slice along the first direction
  if (nx / 2 > nx_local_offset && nx / 2 <= nx_local_offset + nx_local)
    owner = mype;

  if (dir == 1) {
    ix = nx / 2 - 1; // Selecting the slice in the middle

    if (mype == owner) {
      fp = fopen(buf, "w");

      for (iy = 0; iy < ny; ++iy) { // Print on file along the first dimension
        for (iz = 0; iz < nz; ++iz) {
          index = index_f(ix - nx_local_offset, iy, iz, ny, nz);
          fprintf(fp, "%14.6f", data[index]);
        }

        fprintf(fp, "\n");
      }

      fclose(fp);
    }
  } else if (dir == 2 || dir == 3) {
    // In these 2 cases the slices are scattered
    iy = ny / 2 - 1;
    iz = nz / 2 - 1;
    int nk = (dir == 2) ? nz : ny;
    int ik;

    sizes = (int *)malloc(npes * sizeof(int)); // Contains each local dimension
    displ = (int *)calloc(npes, sizeof(int));  // Contains slice along direction
    buffer = (double *)malloc(nx * nk * sizeof(double));
    buffer1d = (double *)malloc(nx /* BUG: nk */ * sizeof(double));

    local_buffer = (double *)malloc(nx_local * sizeof(double));

    MPI_Gather(&nx_local, 1, MPI_INT, sizes, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (mype == 0) { // Only the root computes the displacements
      for (i = 1; i < npes; ++i)
        displ[i] = sizes[i - 1] + displ[i - 1];
    }

    if (dir == 2) { // Select one direction
      for (iz = 0; iz < nz; ++iz) {
        for (ix = 0; ix < nx_local; ++ix) {
          index = index_f(ix, iy, iz, ny, nz);
          local_buffer[ix] = data[index];
        }

        // Gathering each vertical colum
        MPI_Gatherv(local_buffer, nx_local, MPI_DOUBLE, buffer1d, sizes, displ, MPI_DOUBLE, 0,
                    MPI_COMM_WORLD);

        for (ix = 0; ix < nx; ++ix) {
          buffer[ix * nz + iz] = buffer1d[ix]; // Copying each vertical column
        }
      }
    } else { // Or the other
      for (iy = 0; iy < ny; ++iy) {
        for (ix = 0; ix < nx_local; ++ix) {
          index = index_f(ix, iy, iz, ny, nz);
          local_buffer[ix] = data[index];
        }

        // Gathering each vertical colum
        MPI_Gatherv(local_buffer, nx_local, MPI_DOUBLE, buffer1d, sizes, displ, MPI_DOUBLE, 0,
                    MPI_COMM_WORLD);

        for (ix = 0; ix < nx; ++ix) {
          buffer[ix * ny + iy] = buffer1d[ix]; // Copying each vertical column
        }
      }
    }

    if (mype == 0) {
      fp = fopen(buf, "w");

      for (ix = 0; ix < nx; ++ix) {
        for (ik = 0; ik < nk; ++ik) {
          fprintf(fp, "%14.6f", buffer[ix * nk + ik]);
        }

        fprintf(fp, "\n");
      }

      fclose(fp);
    }

    free(sizes);
    free(displ);
    free(buffer);
    free(buffer1d);
    free(local_buffer);
  } else
    fprintf(stderr, " Wrong value for argument 7 in plot_data_2d \n");
}
