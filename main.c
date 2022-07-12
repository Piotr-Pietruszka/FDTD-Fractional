#include "lib_fdtd_fractional.h"

int main()
{

    double alpha_array[14] = {0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 0.995, 0.999, 0.9995, 0.9999};
    for(int i = 0; i < 14; i++)
    {
    // domain constants
    double dz = 0.02e-6;
#ifdef FRACTIONAL_SIM
    double alpha = alpha_array[i];
    double dt_analytical = pow(2.0, 1.0-1.0/alpha) * pow(sqrt(EPS_0*MU_0) * dz, 1.0/alpha);
    double dt = 0.999*dt_analytical; //3.0757e-17; // 2.3068e-17
#else
    double alpha = 1.0;
    double dt = 0.999*dz/C_CONST; 
#endif

    double Lz = 40.0e-6;
    // double Lz = 190.0e-6;
    int Nz = (int) (Lz/dz);
    double T = 6e-14;
    // double T = 20e-14;
    int Nt =  (int) (T/dt);
    
    printf("alpha= %f\n", alpha);
    printf("Nz= %d\n", Nz);
    printf("Nt= %d\n", Nt);
    printf("dz= %e\n", dz);
    printf("dt= %e\n", dt);

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
    int k_source = (int) (20e-6/dz);//(int) (3.03e-6/dz);
    printf("k_source= %d\n", k_source);

    for (int t = 0; t < Nt-1; t++)
    {
        double delay = 7.957747154594768e-15 - dt/2.0;
        // Ex_source[t] = 1/0.3135*1.0/1.444*dt/2.3281e-17*cos((t*dt-delay)*3.707079331235956e+15) * exp( -pow((t*dt-delay) / (1.989436788648692e-15), 2.0) ); // modulated gaussian
        Ex_source[t] = 1.9833*pow(dt, alpha)/6.664611e-017*cos((t*dt-delay)*3.707079331235956e+15) * exp( -pow((t*dt-delay) / (1.989436788648692e-15), 2.0) ); // modulated gaussian

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
    saveFieldToBinary(filename, Ex_source, 1, Nt, dz, dt, alpha);

    double start_time, end_time;
    start_time = omp_get_wtime();
    
    // checkStability();
    simulation(dz, Nz, dt, Nt, alpha, Ex, Hy, Ex_source, k_source, 1);

    end_time = omp_get_wtime();
    double sim_time = end_time - start_time;
    printf("simulation time: %lf s\n", sim_time);

    sprintf(filename, ".\\results\\results.bin");
    saveSimParamsToBinary(filename, dz, Lz, Nz,
                        dt, T, Nt, alpha, sim_time);

    // free allocated memory
    free(Ex);
    free(Hy);
    free(Ex_source);
    }

}