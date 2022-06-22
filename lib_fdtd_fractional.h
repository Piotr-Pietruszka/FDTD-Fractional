
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

// #define _OPENMP
#include <omp.h>

#define MU_0 1.2566371e-6
#define EPS_0 8.85418781762e-12
#define C_CONST 299792458.0

void HyUpdate(const double dz, const int Nz, const int k_bound, const double dt, const int Nt,
                const double* Ex, double* Hy, const int t, const double alpha, 
                const double* GL_coeff_arr);

void ExUpdate(const double dz, const int Nz, const int k_bound, const double dt, const int Nt,
                double* Ex, const double* Hy, const int t, const double alpha, 
                const double* GL_coeff_arr);

void HyClassicUpdate(const double dz, const int Nz, const int k_bound, const double dt, const int Nt,
                     const double* Ex, double* Hy, const int t);

void ExClassicUpdate(const double dz, const int Nz, const int k_bound, const double dt, const int Nt,
                     double* Ex, const double* Hy, const int t);

void simulation(const double dz, const int Nz, const double dt, const int Nt,
                const double alpha,
                double* Ex, double* Hy,
                double* Ex_source, const int k_source, 
                int save_result);

void saveFieldToBinary(const char *filename,
							const double *data,
							const unsigned int Nz,
							const unsigned int Nt,
                            const double dz,
                            const double dt);

double binomialCoeff(const double alpha, const int k);

double fracGLCoeff(const double w, const double alpha, const int n);

void checkStability();

double findMaxAbsValue(const int Nz, const int Nt, double* Ex, double* Hy);


#define ALPHA_ST_NUM 14
#define DT_ST_NUM 16