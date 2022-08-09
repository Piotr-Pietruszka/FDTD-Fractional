
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <omp.h>

#define MU_0 1.2566371e-6
#define EPS_0 8.85418781762e-12
#define C_CONST 299792458.0

#define ALPHA_ST_NUM 14
#define DT_ST_NUM 14


#define FRACTIONAL_SIM
#define MUR_CONDITION
#define TIME_ROW_WISE // faster, default option

// #define STABILITY_CHECK
#define ADD_SOURCE

#ifdef FRACTIONAL_SIM
    #define OPEN_MP_SPACE
#else
    #define OPEN_MP_SPACE
#endif


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

double simulation(const double dz, const int Nz, const double dt, const int Nt,
                  const double alpha,
                  double* Ex, double* Hy,
                  double* Ex_source, const int k_source, 
                  int save_result);

void saveSimParamsToTxt(const char *filename,
                           const double dz, const double Lz, const unsigned int Nz,
                           const double dt, const double T, const unsigned int Nt,
                           const double alpha, const double sim_time);

void saveSimParamsToBinary(const char *filename,
                           const double dz, const double Lz, const unsigned int Nz,
                           const double dt, const double T, const unsigned int Nt,
                           const double alpha, const double sim_time);

void saveFieldToBinary(const char *filename,
							const double *data,
							const unsigned int Nz,
							const unsigned int Nt,
                            const double dz,
                            const double dt,
                            const double alpha);

double binomialCoeff(const double alpha, const int k);

double fracGLCoeff(const double w, const double alpha, const int n);

void checkStability();

double findMaxAbsValue(const int Nz, const int Nt, double* Ex);

double findMaxAbsValueLastTimeStep(const int Nz, const int Nt, double* Ex);
