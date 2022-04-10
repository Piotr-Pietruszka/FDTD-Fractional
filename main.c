#include "lib_fdtd_fractional.h"

int main()
{

    // domain constants
    double dz = 0.1e-6;

    double dt = 0.1*0.99*dz/C_CONST;
    // double dt = 1e-18;

    double Lz = 250.0e-6;
    int Nz = (int) (Lz/dz);
    double T = 20e-14;
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

    simulation(dz, Nz, dt, Nt, alpha, Ex, Hy, Ex_inc, Hy_inc);

}