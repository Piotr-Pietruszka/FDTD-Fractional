#include "lib_fdtd_fractional.h"

int main()
{

    // domain constants
    double dz = 0.07e-6;

    double dt = 0.1*0.99*dz/C_CONST;
    // double dt = 1e-18;

    double Lz = 140.0e-6;
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

    double* Ex_source = calloc(Nt, sizeof(double));
    

    // Source
    for (int t = 0; t < Nt-1; t++)
    {
        Ex_source[t] = sin(t*dt*2*3.14/0.175e-14) * exp( -pow((t*dt-0.75e-14) / (0.2e-14), 2.0) );
    }
    char filename[128];
    sprintf(filename, ".\\results\\source.bin");
    saveFieldToBinary(filename, Ex_source, 1, Nt);

    double start_time, end_time;
    start_time = omp_get_wtime();
    // return 1;
    simulation(dz, Nz, dt, Nt, alpha, Ex, Hy, Ex_source);

    end_time = omp_get_wtime();
    double sim_time = end_time - start_time;
    printf("simulation time: %lf s\n", sim_time);

}