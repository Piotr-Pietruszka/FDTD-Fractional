#include "lib_fdtd_fractional.h"

enum SourceType {MODULATED_GAUSSIAN, TRIANGLE, RECTANGLE, GAUSSIAN, SINUS};

int main()
{   

    // parameters to choose: alpha, T, Lz, source type
    // domain constants
    // double dz = 0.005e-6;
    double dz = 0.01e-6;
    // double dz = 0.005e-6;

#ifdef FRACTIONAL_SIM
    double alpha = 0.99;
    double dt_analytical = pow(2.0, 1.0-1.0/alpha) * pow(sqrt(EPS_0*MU_0) * dz, 1.0/alpha);
    double dt = 0.999*dt_analytical; 
    double dt_cl_ideal = dz/C_CONST;
#else
    double alpha = 1.0;
    double dt = 0.999*dz/C_CONST; 
#endif

    // double Lz = 140.0e-6;
    double Lz = 75.0e-6;
    double T = 20e-14;
    // double T = 20e-14;
    char source_char = 'm';

    unsigned int Nz = (int) (Lz/dz);
    unsigned int Nt =  (int) (T/dt);

    printf("alpha= %f\n", alpha);
    printf("Lz = %e\n", Lz);
    printf("T = %e\n", T); 
    printf("dz= %e\n", dz);
    printf("dt= %e\n", dt);
    printf("Nz= %d\n", Nz);
    printf("Nt= %d\n", Nt);

    // field arrays - for whole domain and simulation time
    double* Ex = calloc(Nz*Nt, sizeof(double));
    double* Hy = calloc(Nz*Nt, sizeof(double));
    // source in time
    double* Ex_inc = calloc(Nt, sizeof(double));
    double* Hy_inc = calloc(Nt, sizeof(double));

    if(Ex==NULL || Hy == NULL || Ex_inc == NULL || Hy_inc == NULL) 
    {
        printf("Error - couldn't allocate memory for arrays\n"); 
        exit(1); 
    } 

    // Source
    int k_source = (int) (10e-6/dz);
    printf("k_source= %d\n", k_source);

    double max_fr = 750e12;
    double min_fr = 430e12;
    double central_fr = (max_fr+min_fr) / 2.0;
    double w_central = 2* M_PI * central_fr;

    double tau = 2.0 / M_PI / (max_fr-min_fr);
    double delay = 4*tau+dt/2;
    for (int t = 0; t < Nt-1; t++)
    {
        Hy_inc[t] = -sqrt(EPS_0/MU_0) * cos(((t+0.5+0.5*dt_cl_ideal/dt)*dt-delay)*w_central) * exp( -pow(((t+0.5+0.5*dt_cl_ideal/dt)*dt-delay) / tau, 2.0) );
        Ex_inc[t] = 1.0 * cos((t*dt-delay)*w_central) * exp( -pow((t*dt-delay) / tau, 2.0) );
    }

    char filename[128];
    sprintf(filename, ".\\results\\source.bin");
    saveFieldToBinary(filename, Ex_inc, 1, Nt, dz, dt, alpha, Nz+10);

    int k_bound = k_source + (int) (20e-6/dz); // material boundary
    double sim_time = simulation(dz, Nz, dt, Nt, alpha, Ex, Hy, Ex_inc, Hy_inc, k_source, k_bound, 1);
    printf("simulation time: %lf s\n", sim_time);

    sprintf(filename, ".\\results\\results.txt");
    saveSimParamsToTxt(filename, dz, Lz, Nz,
                        dt, T, Nt, alpha, sim_time);

    // free allocated memory
    free(Ex);
    free(Hy);
    free(Ex_inc);
    free(Hy_inc);
}