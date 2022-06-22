#include "lib_fdtd_fractional.h"

int main()
{

    // domain constants
    // double dz = 0.03e-6;
    double dz = 0.02e-6;

    // double dt = 0.99*dz/C_CONST;
    // double dt = 0.25*0.99*dz/C_CONST;
    // double dt = 0.6*0.99*dz/C_CONST;

    double dt = 2.0429e-17; //3.0757e-17; // 2.3068e-17

    // double dt = 0.2*0.99*dz/C_CONST; // 0.99
    // double dt = 1e-18;

    double Lz = 150.0e-6;
    // double Lz = 190.0e-6;
    int Nz = (int) (Lz/dz);
    double T = 25e-14;
    // double T = 20e-14;
    // double T = 1.3307e-14;
    int Nt =  (int) (T/dt);
    
    double alpha = 0.99;
    // double alpha = 1.0;
    
    printf("alpha= %f\n", alpha);
    printf("Nz= %d\n", Nz);
    printf("Nt= %d\n", Nt);
    printf("dz= %e\n", dz);
    printf("dt= %e\n", dt);
    printf("k_source= %d\n", (int) (3.03e-6/dz));

    // TEMP - init conditions
    // Nt =  20;

    // field arrays - for whole domain and simulation time
    double* Ex = calloc(Nz*Nt, sizeof(double));
    double* Hy = calloc(Nz*Nt, sizeof(double));
    // source in time
    double* Ex_source = calloc(Nt, sizeof(double));

    if(Ex==NULL || Hy == NULL || Ex_source == NULL) 
    {
        printf("Error - couldn't allocate memory for arrays\n"); 
        exit(1); 
    } 
    
    // TEMP - init conditions
    // Ex[0+ 501*Nt] = 1.0;

    // Source
    int k_source = (int) (3.03e-6/dz);
    for (int t = 0; t < Nt-1; t++)
    {
        double delay = 7.957747154594768e-15 - dt/2.0;
        Ex_source[t] = 1/0.3135*1.0/1.444*dt/2.3281e-17*cos((t*dt-delay)*3.707079331235956e+15) * exp( -pow((t*dt-delay) / (1.989436788648692e-15), 2.0) ); // modulated gaussian

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
    
    // checkStability();
    simulation(dz, Nz, dt, Nt, alpha, Ex, Hy, Ex_source, k_source, 1);

    end_time = omp_get_wtime();
    double sim_time = end_time - start_time;
    printf("simulation time: %lf s\n", sim_time);

    // free allocated memory
    free(Ex);
    free(Hy);
    free(Ex_source);

}