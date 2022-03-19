
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// #define _OPENMP
#include <omp.h>

#define MU_0 1.2566371e-6
#define EPS_0 8.85418781762e-12
#define C_CONST 299792458.0

void HyUpdate(const double dz, const int Nz, const double dt, const int Nt,
                const double* Ex, double* Hy, const int t, const double alpha);

void ExUpdate(const double dz, const int Nz, const double dt, const int Nt,
                double* Ex, const double* Hy, const int t, const double alpha);

void simulation(const double dz, const int Nz, const double dt, const int Nt,
                const double alpha,
                double* Ex, double* Hy);

void saveFieldToBinary(const char *filename,
							const double *data,
							const unsigned int Nz,
							const unsigned int Nt);

double binomialCoeff(const double alpha, const int k);

double wCoeff(const double w, const double alpha, const int n);