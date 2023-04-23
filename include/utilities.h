// Created by G.P. Brandino, I. Girotto, R. Gebauer
// Modified by W. Nadalin
// Last revision: April 2023

#ifndef _FFTW_UTLITIES_
#define _FFTW_UTLITIES_
#include <complex.h>
#include <math.h>
#include <mpi.h>
#include <fftw3.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#ifdef _FFTW3_MPI
#include <fftw3-mpi.h>
#endif

#define pi 3.14159265358979323846

#define MAX(x, y) (((x) > (y)) ? (x) : (y))

typedef struct {
#ifdef _FFTW3_MPI
  fftw_plan fw_plan;
  fftw_plan bw_plan;
#else
  fftw_plan fw_plan_1d;
  fftw_plan fw_plan_2d;

  fftw_plan bw_plan_1d;
  fftw_plan bw_plan_2d;

  fftw_complex *buffer;

  int *sendcounts;
  int *recvcounts;
  int *rdispls;
  int *sdispls;

  ptrdiff_t nx;
  ptrdiff_t ny;
  ptrdiff_t nz;

  ptrdiff_t local_size_slice;
  ptrdiff_t local_ny;

  MPI_Datatype *sendtypes;
  MPI_Datatype *recvtypes;

  MPI_Datatype blocks;
#endif
  fftw_complex *data;

  ptrdiff_t global_size_grid;
  ptrdiff_t local_size_grid;
  ptrdiff_t local_nx;
  ptrdiff_t local_nx_offset;

  MPI_Comm mpi_comm;
} fftw_dist_handler;

double seconds();

int index_f(int, int, int, int, int);

void plot_data_2d(char *, int, int, int, int, int, int, double *);
void print_info(double *, int, int, int, int, int, double, double, double, int, int, double);

void init_fftw(fftw_dist_handler *, int, int, int, MPI_Comm);
void close_fftw(fftw_dist_handler *);

void derivative(fftw_dist_handler *, int, int, int, double, double, double, int, double *,
                double *);
void fft_3d(fftw_dist_handler *, double *, fftw_complex *, bool);

#endif
