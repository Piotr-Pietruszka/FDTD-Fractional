#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <omp.h>

/* ------------------ Physical constants ------------------*/
#define MU_0 1.2566371e-6
#define EPS_0 8.85418781762e-12
#define C_CONST 299792458.0

/* ------------------ Simulation flags and options ------------------*/
#define OPEN_MP_SPACE
#define MUR_CONDITION
#define TIME_ROW_WISE // faster, default option for indexation of field arrays

#define FRACTIONAL_SIMULATION 0
#define DIFFERENT_MATERIALS 1
#define CLASSICAL_SIMULATION 2
#define SIMULATION_TYPE DIFFERENT_MATERIALS

enum SourceType {MODULATED_GAUSSIAN, TRIANGLE, RECTANGLE, GAUSSIAN};

#define CHUNK_SIZE_BYTES 2147483648 // 2 GB (2147483648 bytes); used in saving file to bin

/* ------------------ Function declarations ------------------*/
void initializeSource(const double dz, const double dt, const int Nt, const double alpha,
                      double* Ex_source, enum SourceType source_type);

void initializeSourceTfsf(const double dt, const int Nt, const double dt_cl_ideal,
                          double* Ex_source,  double* Hy_source, enum SourceType source_type);

void HyUpdate(const double dz, const int Nz, const double dt, const int Nt,
              const double* Ex, double* Hy, const int t, const double alpha, 
              const double* GL_coeff_arr);

void ExUpdate(const double dz, const int Nz, const double dt, const int Nt,
              double* Ex, const double* Hy, const int t, const double alpha, 
              const double* GL_coeff_arr);

void HyUpdateDiffMaterials(const double dz, const int Nz, const int k_bound, 
                           const double dt, const int Nt,
                           const double* Ex, double* Hy, const int t, const double alpha,
                           const double* GL_coeff_arr);

void ExUpdateDifferentMaterials(const double dz, const int Nz, const int k_bound, 
                                const double dt, const int Nt,
                                double* Ex, const double* Hy, const int t, const double alpha,
                                const double* GL_coeff_arr);

void HyUpdateClassical(const double dz, const int Nz, const double dt, const int Nt,
                       const double* Ex, double* Hy, const int t);

void ExUpdateClassical(const double dz, const int Nz, const double dt, const int Nt,
                       double* Ex, const double* Hy, const int t);

double simulation(const double dz, const int Nz, const double dt, const int Nt,
                  const double alpha,
                  double* Ex, double* Hy,
                  double* Ex_source, double* Hy_source, const int k_source, 
                  const int k_bound, int save_result);

void saveSimParamsToTxt(const char *filename,
                           const double dz, const double Lz, const unsigned int Nz,
                           const double dt, const double T, const unsigned int Nt,
                           const double alpha, const double sim_time);

void saveFieldToBinary(const char *filename,
							const double *data,
							const unsigned int Nz,
							const unsigned int Nt,
                            const double dz,
                            const double dt,
                            const double alpha, 
                            const int k_bound);

double binomialCoeff(const double alpha, const int k);

double fracGLCoeff(const double w, const double alpha, const int n);
