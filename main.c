#include "lib_fdtd_fractional.h"

int main()
{

    // domain constants
    double dz = 0.1e-6;

    double dt = 0.1*0.99*dz/C_CONST;
    // double dt = 1e-18;

    double Lz = 150.0e-6;
    int Nz = (int) (Lz/dz);
    double T = 25e-14;
    int Nt =  (int) (T/dt);
    double alpha = 0.99;

    printf("Nz: %d\n", Nz);
    printf("Nt: %d\n", Nt);
    printf("dz: %e\n", dz);
    printf("dt: %e\n", dt);

    // field arrays - for whole domain and simulation time
    double* Ex = calloc(Nz*Nt, sizeof(double));
    double* Hy = calloc(Nz*Nt, sizeof(double));

    double* Ex_inc = calloc(Nt, sizeof(double));
    double* Hy_inc = calloc(Nt, sizeof(double));

    clock_t start_time, end_time;
    start_time = clock();

    simulation(dz, Nz, dt, Nt, alpha, Ex, Hy, Ex_inc, Hy_inc);

    end_time = clock();
    double sim_time = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;
    printf("simulation time: %lf\n", sim_time);

}