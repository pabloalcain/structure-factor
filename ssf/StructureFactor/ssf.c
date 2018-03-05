#include "ssf.h"

void ssf(double *x, int natoms, double size, int npoints, int naver, int nrep,
         double *k, double *sk){
  /* Calculates ssf with in 3D for a cubic box and returns on ssf. x is
     an array that has the info in XYZ XYZ XYZ...fashion. The rdf is
     calculated for the whole box up to sqrt(3)/2.


     on sk, returns a 2D array with this format:

     column 0: k
     column k: s(k) in the pair labeled by k
  */
  int hist = nrep * nrep * nrep * natoms;
  
  for (int i = 0; i < 2 * npoints; i++)
    sk[i] = 0;

  double *qx = (double *) malloc(naver*sizeof(double));
  double *qy = (double *) malloc(naver*sizeof(double));
  double *qz = (double *) malloc(naver*sizeof(double));
  double *qw = (double *) malloc(naver*sizeof(double));
  ld_by_order(naver, qx, qy, qz, qw);
#ifdef _OPENMP 
#pragma omp parallel
#endif
  {
    int tid = 0;
    int nthreads = 1;
#ifdef _OPENMP
    tid = omp_get_thread_num();
    nthreads = omp_get_num_threads();
#endif
    for (int it = 0; it < npoints; it += nthreads) {
      int ii = it + tid;
      if (ii >= npoints) break;
      double s_cont = 0.0;
      double ki = k[ii];
      sk[2 * ii] = ki;
      /* average in the sphere */
      for (int j = 0; j < naver; j++) {
        double q1[3];
        q1[0] = ki * qx[j];
        q1[1] = ki * qy[j];
        q1[2] = ki * qz[j];

        /* sum over all atom positions */
        double cell_real = 0.0, cell_imag = 0.0;

        for (int i = 0; i < natoms; i++) {
          double accum;
          accum = q1[0] * x[3*i + 0];
          accum += q1[1] * x[3*i + 1];
          accum += q1[2] * x[3*i + 2];
          double this_real, this_imag;
          this_real = cos(accum);
          this_imag = sin(accum);
          cell_real += this_real;
          cell_imag += this_imag;
        }

        /* sum over all cell images */
        double pbc_real, pbc_imag;
        pbc_real = 0.0;
        pbc_imag = 0.0;
        for (int i = 0; i < nrep; i++) {
          for (int l = 0; l < nrep; l++) {
            for (int m = 0; m < nrep; m++) {
              double accum;
              accum = q1[0] * size * i;
              accum += q1[1] * size * l;
              accum += q1[2] * size * m;
              pbc_real += cos(accum);
              pbc_imag += sin(accum);
            }
          }
        }
        double contr = cell_real * cell_real + cell_imag * cell_imag;
        contr *= pbc_real * pbc_real + pbc_imag * pbc_imag;
        s_cont += contr * qw[j];
      }
      sk[2 * ii + 1] = s_cont;
    }
  }

  for (int i = 0; i < npoints; i++) {
    sk[2 * i + 1] /= hist;
  }
  free(qx);
  free(qy);
  free(qz);
  free(qw);
  return;
}
