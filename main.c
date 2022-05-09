#include "lib_fdtd_fractional.h"

int main()
{

    // domain constants
    double dz = 0.03e-6;

    double dt = 0.2*0.99*dz/C_CONST;
    // double dt = 0.2*0.99*dz/C_CONST; // 0.99
    // double dt = 1e-18;

    double Lz = 170.0e-6;
    int Nz = (int) (Lz/dz);
    // double T = 25e-14;
    double T = 5e-14;
    int Nt =  (int) (T/dt);
    double alpha = 0.98;
    // double alpha = 1.0;

    printf("Nz= %d\n", Nz);
    printf("Nt= %d\n", Nt);
    printf("dz= %e\n", dz);
    printf("dt= %e\n", dt);

    // field arrays - for whole domain and simulation time
    double* Ex = calloc(Nz*Nt, sizeof(double));
    double* Hy = calloc(Nz*Nt, sizeof(double));

    double* Ex_source = calloc(Nt, sizeof(double));
    

    // Source
    for (int t = 0; t < Nt-1; t++)
    {
        // Ex_source[t] = sin(t*dt*2*3.14/0.175e-14) * exp( -pow((t*dt-0.75e-14) / (0.2e-14), 2.0) ); // modulated gaussian
        Ex_source[t] = cos((t*dt-0.75e-14)*2*3.14/0.175e-14) * exp( -pow((t*dt-0.75e-14) / (0.2e-14), 2.0) ); // modulated gaussian
        // Ex_source[t] = (t*dt > 0.3e-14 && t*dt < 0.7e-14) ? 1.0 : 0.0; // rectangle - doesn't work

        // Triangle
        // if(t*dt > 0.4e-14 && t*dt < 0.5e-14)
        // {
        //     Ex_source[t] = (t*dt - 0.4e-14) * 1 / 0.1e-14; // linear function growing to reach 1 at 0.5
        // }
        // else if(t*dt > 0.5e-14 && t*dt < 0.6e-14)
        // {
        //     Ex_source[t] = 1.0 - (t*dt - 0.5e-14) * 1 / 0.1e-14; // linear function decreasing to reach 0 at 0.6
        // }
        
    }
    char filename[128];
    sprintf(filename, ".\\results\\source.bin");
    saveFieldToBinary(filename, Ex_source, 1, Nt, dz, dt);

    double start_time, end_time;
    start_time = omp_get_wtime();
    // return 1;
    simulation(dz, Nz, dt, Nt, alpha, Ex, Hy, Ex_source);

    end_time = omp_get_wtime();
    double sim_time = end_time - start_time;
    printf("simulation time: %lf s\n", sim_time);

}