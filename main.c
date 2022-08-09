#include "lib_fdtd_fractional.h"

enum SourceType {MODULATED_GAUSSIAN, TRIANGLE, RECTANGLE};

int main()
{   


    // domain constants
    double dz = 0.01e-6;
#ifdef FRACTIONAL_SIM
    double alpha = 0.98;
    double dt_analytical = pow(2.0, 1.0-1.0/alpha) * pow(sqrt(EPS_0*MU_0) * dz, 1.0/alpha);
    double dt = 0.9999*dt_analytical; //3.0757e-17; // 2.3068e-17
#else
    double alpha = 1.0;
    double dt = 0.999*dz/C_CONST; 
#endif

    double Lz = 50.0e-6;
    // double Lz = 190.0e-6;
    unsigned int Nz = (int) (Lz/dz);
    double T = 5e-14;
    // double T = 20e-14;
    unsigned int Nt =  (int) (T/dt);
    
    printf("alpha= %f\n", alpha);
    printf("Nz= %d\n", Nz);
    printf("Nt= %d\n", Nt);
    printf("dz= %e\n", dz);
    printf("dt= %e\n", dt);
    printf("Lz = %e\n", Lz);
    printf("T = %e\n", T); 

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

    enum SourceType source_type = MODULATED_GAUSSIAN;
    // Source
    int k_source = (int) (0.1e-6/dz);
    printf("k_source= %d\n", k_source);


    // First and last timesteps of source signal
    int first_t_source = -1;
    int last_t_source = -1;
   
    for (int t = 0; t < Nt-1; t++)
    {
        double delay = 7.957747154594768e-15 - dt/2.0;
        if (source_type == MODULATED_GAUSSIAN)
        {
            Ex_source[t] = 5.9916e+08 * pow(dt, alpha) / dz * cos((t*dt-delay)*3.707079331235956e+15) * exp( -pow((t*dt-delay) / (1.989436788648692e-15), 2.0) ); // modulated gaussian       
        }
        else if (source_type == TRIANGLE)
        {
            if(t*dt > 0.4e-14 && t*dt < 0.5e-14)
            {
                Ex_source[t] = (t*dt - 0.4e-14) * 1 / 0.1e-14; // linear function growing to reach 1 at 0.5
            }
            else if(t*dt > 0.5e-14 && t*dt < 0.6e-14)
            {
                Ex_source[t] = 1.0 - (t*dt - 0.5e-14) * 1 / 0.1e-14; // linear function decreasing to reach 0 at 0.6
            }
        }
        else if(source_type == RECTANGLE)
        {
            if(t*dt > 0.3e-14 && t*dt < 0.7e-14)
            {
                if(first_t_source < 0)
                    first_t_source = t;
                Ex_source[t] = 2.0; // rectangular function
                last_t_source = t;
            }
        }
    }
    if(source_type == RECTANGLE)
    {
        Ex_source[first_t_source-1] = 1.75; Ex_source[first_t_source-2] = 1.5; Ex_source[first_t_source-3] = 1.25; Ex_source[first_t_source-4] = 1; Ex_source[first_t_source-5] = 0.75; Ex_source[first_t_source-6] = 0.5; Ex_source[first_t_source-7] = 0.25;
        Ex_source[last_t_source+1] = 1.5; Ex_source[last_t_source+2] = 1; Ex_source[last_t_source+3] = 0.5;
    }

    char filename[128];
    sprintf(filename, ".\\results\\source.bin");
    saveFieldToBinary(filename, Ex_source, 1, Nt, dz, dt, alpha);

#ifdef STABILITY_CHECK
    checkStability();
#else
    double sim_time = simulation(dz, Nz, dt, Nt, alpha, Ex, Hy, Ex_source, k_source, 1);
    printf("simulation time: %lf s\n", sim_time);

    // sprintf(filename, ".\\results\\results.txt");
    // saveSimParamsToTxt(filename, dz, Lz, Nz,
    //                     dt, T, Nt, alpha, sim_time);
#endif

    // free allocated memory
    free(Ex);
    free(Hy);
    free(Ex_source);
}