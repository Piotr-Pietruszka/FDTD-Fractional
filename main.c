#include "lib_fdtd_fractional.h"

int main()
{

    // domain constants
    // double dz = 0.03e-6;
    double dz = 0.02e-6;

    // double dt = 0.99*dz/C_CONST;
    // double dt = 0.25*0.99*dz/C_CONST;
    // double dt = 0.6*0.99*dz/C_CONST;
    double dt = 3.0857e-17;
    // double dt = 0.2*0.99*dz/C_CONST; // 0.99
    // double dt = 1e-18;

    double Lz = 150.0e-6;
    // double Lz = 190.0e-6;
    int Nz = (int) (Lz/dz);
    double T = 35e-14;
    // double T = 20e-14;
    // double T = 1.3307e-14;
    int Nt =  (int) (T/dt);
    
    // double alpha = 0.97;
    double alpha = 0.98;
    
    printf("alpha= %f\n", alpha);
    printf("Nz= %d\n", Nz);
    printf("Nt= %d\n", Nt);
    printf("dz= %e\n", dz);
    printf("dt= %e\n", dt);
    printf("k_source= %d\n", (int) (3.03e-6/dz));

    // field arrays - for whole domain and simulation time
    double* Ex = calloc(Nz*Nt, sizeof(double));
    double* Hy = calloc(Nz*Nt, sizeof(double));

    double* Ex_source = calloc(Nt, sizeof(double));
    
    // TEMP - init conditions
    // Nt =  20;
    // Ex[0+ 501*Nt] = 1.0;

    // Source
    for (int t = 0; t < Nt-1; t++)
    {
        double delay = 7.957747154594768e-15 - dt/2.0;
        Ex_source[t] = 1/3.2*1/0.4624*1/1.473*0.9983/1.002*dt/2.3281e-17*cos((t*dt-delay)*3.707079331235956e+15) * exp( -pow((t*dt-delay) / (1.989436788648692e-15), 2.0) ); // modulated gaussian

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

    simulation(dz, Nz, dt, Nt, alpha, Ex, Hy, Ex_source);

    end_time = omp_get_wtime();
    double sim_time = end_time - start_time;
    printf("simulation time: %lf s\n", sim_time);

}