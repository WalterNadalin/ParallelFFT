// Created by G.P. Brandino, I. Girotto, R. Gebauer
// Modified by W. Nadalin
// Last revision: April 2023

#include "utilities.h"

double seconds() {
  // Return the second elapsed since Epoch (00:00:00 UTC, January 1, 1970)
  struct timeval tmp;
  double sec;
  gettimeofday(&tmp, (struct timezone *)0);
  sec = tmp.tv_sec + ((double)tmp.tv_usec) / 1000000.0;
  return sec;
}

//  Index linearization is computed following row-major order
int index_f(int ix, int iy, int iz, int ny, int nz) { return nz * (ny * ix + iy) + iz; }

void init_fftw(fftw_dist_handler *fft, int nx, int ny, int nz, MPI_Comm comm) {
  fft->mpi_comm = comm;

#ifdef _FFTW3_MPI
  fftw_mpi_init(); // Initialize a parallel enviroment for the fast Fourier transform

  fft->local_size_grid =
      fftw_mpi_local_size_3d(nx, ny, nz, fft->mpi_comm, &(fft->local_nx), &(fft->local_nx_offset));
#else
  int size, rank;

  fft->nx = nx;
  fft->ny = ny;
  fft->nz = nz;

  //  Allocate a distributed grid for complex FFT using aligned memory
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  if (((nx % size) || (ny % size)) && !rank) {
    fprintf(stdout, "\nnx dimension must be multiple of the number of "
                    "processes. The program will be aborted...\n\n");
    MPI_Abort(comm, 1);
  }

  fft->sendcounts = (int *)malloc(size * sizeof(int));
  fft->recvcounts = (int *)malloc(size * sizeof(int));
  fft->rdispls = (int *)malloc(size * sizeof(int));
  fft->sdispls = (int *)malloc(size * sizeof(int));
  fft->sendtypes = (MPI_Datatype *)malloc(size * sizeof(MPI_Datatype));
  fft->recvtypes = (MPI_Datatype *)malloc(size * sizeof(MPI_Datatype));

  // Useful dimensions
  fft->local_nx = fft->nx / size;
  fft->local_ny = fft->ny / size;
  fft->local_nx_offset = fft->local_nx * rank;
  fft->local_size_grid = fft->local_nx * fft->ny * fft->nz;
  fft->local_size_slice = fft->local_nx * fft->local_ny * fft->nz;

  // Datatype which selects the correct data to send
  MPI_Type_vector(fft->local_nx, fft->local_ny * fft->nz, fft->ny * fft->nz, MPI_C_DOUBLE_COMPLEX,
                  &(fft->blocks));
  MPI_Type_commit(&(fft->blocks));

  // Information fo all to all `w` communication
  fft->rdispls[0] = fft->sdispls[0] = 0;
  fft->sendtypes[0] = fft->blocks;
  fft->recvtypes[0] = MPI_C_DOUBLE_COMPLEX;

  for (int i = 0; i < size - 1; i++) {
    fft->sendcounts[i] = 1;
    fft->recvcounts[i] = fft->local_size_slice;
    fft->rdispls[i + 1] = fft->rdispls[i] + fft->local_size_slice * sizeof(fftw_complex);
    fft->sdispls[i + 1] = fft->sdispls[i] + fft->local_ny * fft->nz * sizeof(fftw_complex);
    fft->sendtypes[i + 1] = fft->blocks;
    fft->recvtypes[i + 1] = MPI_C_DOUBLE_COMPLEX;
  }

  fft->sendcounts[size - 1] = 1;
  fft->recvcounts[size - 1] = fft->local_size_slice;

  fft->buffer = (fftw_complex *)fftw_malloc(fft->local_size_grid * sizeof(fftw_complex));
#endif
  fft->global_size_grid = nx * ny * nz;

  // Allocate data to transform and make plans
  fft->data = (fftw_complex *)fftw_malloc(fft->local_size_grid * sizeof(fftw_complex));
#ifdef _FFTW3_MPI
  fft->fw_plan = fftw_mpi_plan_dft_3d(nx, ny, nz, fft->data, fft->data, fft->mpi_comm, FFTW_FORWARD,
                                      FFTW_ESTIMATE);
  fft->bw_plan = fftw_mpi_plan_dft_3d(nx, ny, nz, fft->data, fft->data, fft->mpi_comm,
                                      FFTW_BACKWARD, FFTW_ESTIMATE);
#else
  const int n_2d[] = {ny, nz}, n_1d[] = {nx};

  fft->fw_plan_2d =
      fftw_plan_many_dft(2, n_2d, fft->local_nx, fft->data, n_2d, 1, n_2d[0] * n_2d[1], fft->data,
                         n_2d, 1, n_2d[0] * n_2d[1], FFTW_FORWARD, FFTW_ESTIMATE);

  fft->fw_plan_1d =
      fftw_plan_many_dft(1, n_1d, nz * ny / size, fft->buffer, n_1d, nz * ny / size, 1, fft->buffer,
                         n_1d, nz * ny / size, 1, FFTW_FORWARD, FFTW_ESTIMATE);

  fft->bw_plan_2d =
      fftw_plan_many_dft(2, n_2d, fft->local_nx, fft->data, n_2d, 1, n_2d[0] * n_2d[1], fft->data,
                         n_2d, 1, n_2d[0] * n_2d[1], FFTW_BACKWARD, FFTW_ESTIMATE);

  fft->bw_plan_1d =
      fftw_plan_many_dft(1, n_1d, nz * ny / size, fft->buffer, n_1d, nz * ny / size, 1, fft->buffer,
                         n_1d, nz * ny / size, 1, FFTW_BACKWARD, FFTW_ESTIMATE);
#endif
}

void close_fftw(fftw_dist_handler *fft) {
#ifdef _FFTW3_MPI
  fftw_destroy_plan(fft->bw_plan);
  fftw_destroy_plan(fft->fw_plan);
#else
  fftw_destroy_plan(fft->bw_plan_1d);
  fftw_destroy_plan(fft->bw_plan_2d);

  fftw_destroy_plan(fft->fw_plan_1d);
  fftw_destroy_plan(fft->fw_plan_2d);

  fftw_free(fft->buffer);

  free(fft->sendcounts);
  free(fft->recvcounts);
  free(fft->rdispls);
  free(fft->sdispls);
  free(fft->sendtypes);
  free(fft->recvtypes);

  MPI_Type_free(&(fft->blocks));
#endif
  fftw_free(fft->data);
}

// This subroutine uses fftw to calculate 3-dimensional discrete FFTs
// The data in direct space is assumed to be real-valued
// The data in reciprocal space is complex direct_to_reciprocal indicates in which direction the FFT
// is to be calculated Note that for real data in direct space (like here), we have F(N-j) =
// conj(F(j)) where F is the array in reciprocal space Here, we do not make use of this property
// Also, we do not use the special (time-saving) routines of FFTW which allow one to save time and
// memory for such real-to-complex transforms
// f (array in direct space): f(l) = 1/N \sum_{k=0}^{N-1} exp(+ 2 \pi I k*l/N) F(k)
// F (array in reciprocal space): F(k) = \sum_{l=0}^{N-1} exp(- 2 \pi I k*l/N) f(l)

void fft_3d(fftw_dist_handler *fft, double *data_direct, fftw_complex *data_rec,
            bool direct_to_reciprocal) {
  double fac;
  int i;
  const int local_size_grid = fft->local_size_grid;
  fftw_complex *data = fft->data;

#ifndef _FFTW3_MPI
  const int *sendcounts = fft->sendcounts, *recvcounts = fft->recvcounts, *rdispls = fft->rdispls,
            *sdispls = fft->sdispls;
  MPI_Datatype *sendtypes = fft->sendtypes, *recvtypes = fft->recvtypes;
  fftw_complex *buffer = fft->buffer;
#endif

  if (direct_to_reciprocal) { // Direct transform
    for (i = 0; i < local_size_grid; i++)
      data[i] = data_direct[i] + 0.0 * I;

#ifdef _FFTW3_MPI
    fftw_execute(fft->fw_plan);
#else
    fftw_execute(fft->fw_plan_2d); // 2D transform along planes

    MPI_Alltoallw(data, sendcounts, sdispls, sendtypes, buffer, recvcounts, rdispls, recvtypes,
                  MPI_COMM_WORLD);

    fftw_execute(fft->fw_plan_1d); // 1D transform along columns

    MPI_Alltoallw(buffer, recvcounts, rdispls, recvtypes, data, sendcounts, sdispls, sendtypes,
                  MPI_COMM_WORLD);
#endif

    memcpy(data_rec, data, fft->local_size_grid * sizeof(fftw_complex));
  } else { // Reverse transform
    memcpy(data, data_rec, fft->local_size_grid * sizeof(fftw_complex));

#ifdef _FFTW3_MPI
    fftw_execute(fft->bw_plan);
#else
    fftw_execute(fft->bw_plan_2d); // 2D transform along planes

    MPI_Alltoallw(data, sendcounts, sdispls, sendtypes, buffer, recvcounts, rdispls, recvtypes,
                  MPI_COMM_WORLD);

    fftw_execute(fft->bw_plan_1d); // 1D transform along columns

    MPI_Alltoallw(buffer, recvcounts, rdispls, recvtypes, data, sendcounts, sdispls, sendtypes,
                  MPI_COMM_WORLD);
#endif

    fac = 1.0 / fft->global_size_grid;

    for (i = 0; i < local_size_grid; i++)
      data_direct[i] = creal(data[i]) * fac;
  }
}
