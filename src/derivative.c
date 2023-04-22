// Created by G.P. Brandino, I. Girotto, R. Gebauer
// Modified by W. Nadalin
// Last revision: April 2023

#include "utilities.h"

// Calculate the derivative in direction ipol of the array 'data'
void derivative(fftw_dist_handler *fft, int nx, int ny, int nz, double L1, double L2, double L3,
                int ipol, double *data, double *deriv) {
  fftw_complex *aux;
  double G;
  int i, ix, iy, iz, index;

  /// This implementation will perform a full 3D FFT of the data, and then derive
  aux = (fftw_complex *)fftw_malloc(fft->local_size_grid * sizeof(fftw_complex));

  fft_3d(fft, data, aux, true); // First get the FFT of data

  if (ipol == 1) {

    G = 2.0 * pi / L1;
    for (ix = 0; ix < fft->local_nx; ++ix) {
      i = ix + fft->local_nx_offset;

      if (i > nx / 2)
        i = i - nx;
      if (i == nx / 2)
        i = 0;

      for (iy = 0; iy < ny; ++iy) {
        for (iz = 0; iz < nz; ++iz) {
          index = index_f(ix, iy, iz, ny, nz);
          aux[index] *= 0.0 + G * i * I;
        }
      }
    }
  }

  if (ipol == 2) {
    G = 2.0 * pi / L2;

    for (iy = 0; iy < ny; ++iy) {

      i = iy;
      if (i > ny / 2)
        i = i - ny;
      if (i == ny / 2)
        i = 0;

      for (ix = 0; ix < fft->local_nx; ++ix) {
        for (iz = 0; iz < nz; ++iz) {
          index = index_f(ix, iy, iz, ny, nz);
          aux[index] *= 0.0 + G * i * I;
        }
      }
    }
  }

  if (ipol == 3) {
    G = 2.0 * pi / L3;

    for (iz = 0; iz < nz; ++iz) {

      i = iz;
      if (i > nz / 2)
        i = i - nz;
      if (i == nz / 2)
        i = 0;

      for (ix = 0; ix < fft->local_nx; ++ix) {
        for (iy = 0; iy < ny; ++iy) {
          index = index_f(ix, iy, iz, ny, nz);
          aux[index] *= 0.0 + G * i * I;
        }
      }
    }
  }

  fft_3d(fft, deriv, aux, false); // Now go back to real space
  fftw_free(aux);
}
