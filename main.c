#include "lib_fdtd_fractional.h"

#define NO_DIM 1 // number of different dimension size combinations to run simulation with




/**
 * Init parameters to run simulation with. Goal is to test different
 * parameter values to measure speed depending on Nz, Nt, alpha, ...
 * 
 * 
 * @param w previous coefficient
 * @param alpha fractional order of derivative 
 * @param n coefficient number (index)
 * @return new (next) coefficient
 */
void initParamsCombinations(unsigned int dim_arr[NO_DIM][2])
{

    
    int Nz_change = 0;
    if(Nz_change > NO_DIM)
        Nz_change = NO_DIM;
    
    // Different Nz, constant Nt
    for(int i = 3; i < Nz_change; i++)
    {
        dim_arr[i][0] = i*100;
        dim_arr[i][1] = 300;
    }
    
    // Different Nt, constant Nz
    for(int i = Nz_change; i < NO_DIM; i++)
    {
        dim_arr[i][0] = 4000;
        dim_arr[i][1] = (i+1-Nz_change) * 1000; // Nt value
    }
}





#define RUN_SIM
int main()
{   
    // {{Nz, Nt}, {Nz, Nt}, ...}
    
    int a = omp_get_num_procs();
    printf("%d\n", a);

    // dim_arr[_][0] - Nz, dim_arr[_][1] - Nt
    unsigned int dim_arr[NO_DIM][2] = {100, 100};
    // unsigned int dim_arr[NO_DIM][2] = { {500, 700}, {900, 700} , {1500, 700}, 
    //                                     {2000, 700}, {2500, 700}, {3000, 700}, {3500, 700}, {4000, 700}, {5000, 700}, {7000, 700},
    //                                     {9000, 700}, {11000, 700}, {13000, 700},  {15000, 700}, {17000, 700}};
    // printf("%d, %d, %d, %d\n", dim_arr[0][0], dim_arr[0][1], dim_arr[1][0], dim_arr[2][1]);
    
    // initParamsCombinations(dim_arr);

#ifndef RUN_SIM

    char filename[128];
    sprintf(filename, ".\\results\\results.bin");
    FdtdCpuReadPlanefromFile(filename);

#endif


#ifdef RUN_SIM
for(int i = 0; i < NO_DIM; i++)
    {
    printf("\ni= %d\n", i);    
    // domain constants
    double dz = 0.01e-6;
#ifdef FRACTIONAL_SIM
    double alpha = 0.98;
    double dt_analytical = pow(2.0, 1.0-1.0/alpha) * pow(sqrt(EPS_0*MU_0) * dz, 1.0/alpha);
    double dt = 0.999*dt_analytical; //3.0757e-17; // 2.3068e-17
#else
    double alpha = 1.0;
    double dt = 0.999*dz/C_CONST; 
#endif

    unsigned int Nz = dim_arr[i][0];
    unsigned int Nt = dim_arr[i][1];
    double Lz = Nz*dz;
    printf("Lz = %e\n", Lz); // TEMP
    double T = Nt*dt;
    
    printf("alpha= %f\n", alpha);
    printf("Nz= %d\n", Nz);
    printf("Nt= %d\n", Nt);
    printf("dz= %e\n", dz);
    printf("dt= %e\n", dt);

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

    // Source
    int k_source = (int) (0.1e-6/dz);
    printf("k_source= %d\n", k_source);

    for (int t = 0; t < Nt-1; t++)
    {
        double delay = 7.957747154594768e-15 - dt/2.0;
        Ex_source[t] = 5.9916e+08 * pow(dt, alpha) / dz * cos((t*dt-delay)*3.707079331235956e+15) * exp( -pow((t*dt-delay) / (1.989436788648692e-15), 2.0) ); // modulated gaussian       
    }
    char filename[128];
    sprintf(filename, ".\\results\\source.bin");
    saveFieldToBinary(filename, Ex_source, 1, Nt, dz, dt, alpha);

    // checkStability();

    // double sim_time = simulation(dz, Nz, dt, Nt, alpha, Ex, Hy, Ex_source, k_source, 1);
    double sim_time = simulation(dz, Nz, dt, Nt, alpha, Ex, Hy, Ex_source, k_source, 0);
    printf("simulation time: %lf s\n", sim_time);

    sprintf(filename, ".\\results\\results.txt");
    // saveSimParamsToTxt(filename, dz, Lz, Nz,
    //                     dt, T, Nt, alpha, sim_time);

    // free allocated memory
    free(Ex);
    free(Hy);
    free(Ex_source);

    }
#endif

}