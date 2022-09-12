#include "lib_fdtd_fractional.h"


int main()
{   
#if SIMULATION_TYPE == FRACTIONAL_SIMULATION
    printf("Simulation of fractional order material\n");
#elif SIMULATION_TYPE == DIFFERENT_MATERIALS
    printf("Simulation of different media: vacuum and fractional order material\n");
#elif SIMULATION_TYPE == CLASSICAL_SIMULATION
    printf("Classical simulation (vacuum)\n");
#endif

    // parameters to choose: alpha, T, Lz, source type
    // domain constants
    double dz = 0.01e-6;
    double dt_cl_ideal = dz/C_CONST;

#if SIMULATION_TYPE == CLASSICAL_SIMULATION
    double alpha = 1.0;
    double dt = 0.999*dz/C_CONST; 
#else 
    double alpha = 0.98;
    double dt_analytical = pow(2.0, 1.0-1.0/alpha) * pow(sqrt(EPS_0*MU_0) * dz, 1.0/alpha);
    double dt = 0.999*dt_analytical;
#endif

    enum SourceType source_type = MODULATED_GAUSSIAN;
    double Lz = 110.0e-6;
    double T = 18e-14;
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
    double* Ex_source = calloc(Nt, sizeof(double));
    double* Hy_source = calloc(Nt, sizeof(double)); // used only in tfsf

    if(Ex==NULL || Hy == NULL || Ex_source == NULL) 
    {
        printf("Error - couldn't allocate memory for arrays\n"); 
        exit(1); 
    } 

    // Source and boundary between materials
#if SIMULATION_TYPE == DIFFERENT_MATERIALS
    int k_source = (int) (10e-6/dz);
    printf("k_source= %d\n", k_source);
    int k_bound = k_source + (int) (20e-6/dz); // material boundary
    printf("k_bound= %d\n", k_bound);
    initializeSourceTfsf(dt, Nt, dt_cl_ideal, Ex_source,  Hy_source, source_type);
#else
    int k_source = (int) (0.1e-6/dz);
    printf("k_source= %d\n", k_source);
    int k_bound = -5; // domain is homogeneous
    initializeSource(dz, dt, Nt, alpha, Ex_source, source_type);
#endif
    if(k_source > Nz || k_bound > (int) Nz)
    {
        printf("Source or boundary position exceeds domain size. Exiting");
        exit(1);
    }

    char filename[128];
    sprintf(filename, ".\\results\\source.bin");
    saveFieldToBinary(filename, Ex_source, 1, Nt, dz, dt, alpha, k_bound);

    // Run simulation
    double sim_time = simulation(dz, Nz, dt, Nt, alpha, Ex, Hy, Ex_source, Hy_source, k_source, k_bound, 1);
    printf("simulation execution time: %lf s\n", sim_time);

    // free allocated memory
    free(Ex);
    free(Hy);
    free(Ex_source);
}