#include "lib_fdtd_fractional.h"

/**
 * Update Hy field in 1D domain
 * 
 * @param dz spatial step size
 * @param Nz domain size (cells)
 * @param dt time step
 * @param Nt length of simulation (timesteps)
 * @param Ex Ex field array 
 * @param Hy Hy field array
 * @param t current timestep
 * @param alpha fractional order of derivative
 * @param GL_coeff_arr Nt-size array of Gr-Let derivative coefficients, starts at w1
 * @return 
 */
void HyUpdate(const double dz, const int Nz, const int k_bound, const double dt, const int Nt,
                const double* Ex, double* Hy, const int t, const double alpha,
                const double* GL_coeff_arr)
{

#ifdef TIME_ROW_WISE
    int k = 0;
    int n = 0;
#ifdef OPEN_MP_SPACE
// #pragma omp parallel for
#pragma omp parallel for private(k)
// #pragma omp parallel for collapse(2)
#endif
    for(k = 0; k < Nz-1; k++)
    {
        // update based on Hy rotation
        Hy[t+1 + k*Nt] =  -pow(dt, alpha)/MU_0 * (Ex[t + (k+1)*Nt] - Ex[t + k*Nt])/dz;

        // update based on previous Hy values (from GL derivative)
#ifdef OPEN_MP_TIME
#pragma omp parallel for private(n)
#endif
        for(n = 0; n < t; n++)
        {
            Hy[t+1 + k*Nt] -= Hy[t-n + k*Nt] * GL_coeff_arr[n];
        }
    }
#else
    int k = 0;
    int n = 0;
#ifdef OPEN_MP_SPACE
#pragma omp parallel for private(k)
#endif
    for(k = 0; k < Nz-1; k++)
    {
        // update based on Ex rotation
        Hy[(t+1)*Nz + k] =  -pow(dt, alpha)/MU_0 * (Ex[t*Nz + k+1] - Ex[t*Nz + k])/dz;

        // update based on previous Hy values (from GL derivative)
#ifdef OPEN_MP_TIME
#pragma omp parallel for private(n)
#endif
        for(n = 0; n < t; n++)
        {
            Hy[(t+1)*Nz + k] -= Hy[(t-n)*Nz + k] * GL_coeff_arr[n];
        }
    }
#endif
}


/**
 * Update Ex field in 1D domain
 * 
 * @param dz spatial step size
 * @param Nz domain size (cells)
 * @param dt time step
 * @param Nt length of simulation (timesteps)
 * @param Ex Ex field array 
 * @param Hy Hy field array
 * @param t current timestep
 * @param alpha fractional order of derivative
 * @param GL_coeff_arr Nt-size array of Gr-Let derivative coefficients, starts at w1
 * @return 
 */
void ExUpdate(const double dz, const int Nz, const int k_bound, const double dt, const int Nt,
                double* Ex, const double* Hy, const int t, const double alpha,
                const double* GL_coeff_arr)
{

#ifdef TIME_ROW_WISE
    int k = 1;
    int n = 0;
#ifdef MUR_CONDITION    
    #ifdef OPEN_MP_SPACE
        #pragma omp parallel for private(k)
    #endif
    for(k = 1; k < Nz-1; k++)
    {
        // update based on Hy rotation
        Ex[t+1 + k*Nt] = -pow(dt, alpha)/EPS_0* (Hy[t+1 + k*Nt] - Hy[t+1 + (k-1)*Nt]) / dz; // update based on Hy rotation
        
        // update based on previous Ex values (from GL derivative)
        #ifdef OPEN_MP_TIME
            #pragma omp parallel for private(n)
        #endif
        for(n = 0; n < t; n++)
        {
            Ex[t+1 + k*Nt] -= Ex[t-n + k*Nt] * GL_coeff_arr[n];
        }
    }
    // Mur condition - left boundary
    Ex[t+1 + 0*Nt] = (C_CONST*pow(dt, alpha)-dz)/(C_CONST*pow(dt, alpha)+dz) * Ex[t+1 + 1*Nt] +
                     (C_CONST*pow(dt, alpha))/(C_CONST*pow(dt, alpha)+dz) * (Ex[t + 1*Nt] - Ex[t + 0*Nt]);
    #ifdef OPEN_MP_TIME
        #pragma omp parallel for private(n)
    #endif
    for(n = 0; n < t; n++)
    {
        Ex[t+1 + 0*Nt] -=  (dz)/(C_CONST*pow(dt, alpha)+dz) * GL_coeff_arr[n] * Ex[t-n + 1*Nt];
        Ex[t+1 + 0*Nt] -=  (dz)/(C_CONST*pow(dt, alpha)+dz) * GL_coeff_arr[n] * Ex[t-n + 0*Nt];
    }
    // Mur condition - right boundary   
    Ex[t+1 + (Nz-1)*Nt] = (C_CONST*pow(dt, alpha)-dz)/(C_CONST*pow(dt, alpha)+dz) * Ex[t+1 + (Nz-2)*Nt] +
                          (C_CONST*pow(dt, alpha))/(C_CONST*pow(dt, alpha)+dz) * (Ex[t + (Nz-2)*Nt] - Ex[t + (Nz-1)*Nt]);
    #ifdef OPEN_MP_TIME
        #pragma omp parallel for private(n)
    #endif
    for(n = 0; n < t; n++)
    {
        Ex[t+1 + (Nz-1)*Nt] -=  (dz)/(C_CONST*pow(dt, alpha)+dz) * GL_coeff_arr[n] * Ex[t-n + (Nz-1)*Nt];
        Ex[t+1 + (Nz-1)*Nt] -=  (dz)/(C_CONST*pow(dt, alpha)+dz) * GL_coeff_arr[n] * Ex[t-n + (Nz-2)*Nt];
    }
#else 
    #ifdef OPEN_MP_SPACE
        #pragma omp parallel for private(k)
    #endif
    for(k = 1; k < Nz-1; k++)
    {
        // update based on Hy rotation
        Ex[t+1 + k*Nt] = -pow(dt, alpha)/EPS_0* (Hy[t+1 + k*Nt] - Hy[t+1 + (k-1)*Nt]) / dz; // update based on Hy rotation
        
        // update based on previous Ex values (from GL derivative)
        #ifdef OPEN_MP_TIME
            #pragma omp parallel for private(n)
        #endif
        for(n = 0; n < t; n++)
        {
            Ex[t+1 + k*Nt] -= Ex[t-n + k*Nt] * GL_coeff_arr[n];
        }      
    }
#endif

#else

    int k = 1;
    int n = 0;
#ifdef MUR_CONDITION
    #ifdef OPEN_MP_SPACE
        #pragma omp parallel for private(k)
    #endif
    for(k = 1; k < Nz-1; k++)
    {
        // update based on Hy rotation
        Ex[(t+1)*Nz + k] = -pow(dt, alpha)/EPS_0* (Hy[(t+1)*Nz + k] - Hy[(t+1)*Nz + k-1]) / dz; // update based on Hy rotation
        
        // update based on previous Ex values (from GL derivative)
        #ifdef OPEN_MP_TIME
            #pragma omp parallel for private(n)
        #endif
        for(n = 0; n < t; n++)
        {
            Ex[(t+1)*Nz + k] -= Ex[(t-n)*Nz + k] * GL_coeff_arr[n];
        }
    }
    // Mur condition - left boundary
    Ex[t+1 + 0*Nt] = (C_CONST*pow(dt, alpha)-dz)/(C_CONST*pow(dt, alpha)+dz) * Ex[(t+1)*Nz + 1] +
                     (C_CONST*pow(dt, alpha))/(C_CONST*pow(dt, alpha)+dz) * (Ex[(t)*Nz + 1] - Ex[(t)*Nz + 0*Nt]);
    #ifdef OPEN_MP_TIME
        #pragma omp parallel for private(n)
    #endif
    for(n = 0; n < t; n++)
    {
        Ex[(t+1)*Nz + 0] -=  (dz)/(C_CONST*pow(dt, alpha)+dz) * GL_coeff_arr[n] * Ex[(t-n)*Nz + 1];
        Ex[(t+1)*Nz + 0] -=  (dz)/(C_CONST*pow(dt, alpha)+dz) * GL_coeff_arr[n] * Ex[(t-n)*Nz + 0];
    }
    // Mur condition - right boundary   
    Ex[(t+1)*Nz + Nz-1] = (C_CONST*pow(dt, alpha)-dz)/(C_CONST*pow(dt, alpha)+dz) * Ex[(t+1)*Nz + Nz-2] +
                          (C_CONST*pow(dt, alpha))/(C_CONST*pow(dt, alpha)+dz) * (Ex[(t)*Nz + Nz-2] - Ex[(t)*Nz + Nz-1]);
    #ifdef OPEN_MP_TIME
        #pragma omp parallel for private(n)
    #endif
    for(n = 0; n < t; n++)
    {
        Ex[(t+1)*Nz + Nz-1] -=  (dz)/(C_CONST*pow(dt, alpha)+dz) * GL_coeff_arr[n] * Ex[(t-n)*Nz + Nz-1];
        Ex[(t+1)*Nz + Nz-1] -=  (dz)/(C_CONST*pow(dt, alpha)+dz) * GL_coeff_arr[n] * Ex[(t-n)*Nz + Nz-2];
    }   
#else 
    #ifdef OPEN_MP_SPACE
        #pragma omp parallel for private(k)
    #endif
    for(k = 1; k < Nz-1; k++)
    {
        // update based on Hy rotation
        Ex[(t+1)*Nz + k] = -pow(dt, alpha)/EPS_0* (Hy[(t+1)*Nz + k] - Hy[(t+1)*Nz + k-1]) / dz; // update based on Hy rotation
        
        // update based on previous Ex values (from GL derivative)
        #ifdef OPEN_MP_TIME
            #pragma omp parallel for private(n)
        #endif
        for(n = 0; n < t; n++)
        {
            Ex[(t+1)*Nz + k] -= Ex[(t-n)*Nz + k] * GL_coeff_arr[n];
        }      
    }
#endif

#endif

}


/**
 * Classical (non-fractional) Hy field update in 1D domain
 * 
 * @param dz spatial step size
 * @param Nz domain size (cells)
 * @param dt time step
 * @param Nt length of simulation (timesteps)
 * @param Ex Ex field array 
 * @param Hy Hy field array
 * @param t current timestep
 * @return 
 */
void HyClassicUpdate(const double dz, const int Nz, const int k_bound, const double dt, const int Nt,
                     const double* Ex, double* Hy, const int t)
{
    int k = 1;
#ifdef MUR_CONDITION
    #ifdef OPEN_MP_SPACE
        #pragma omp parallel for private(k)
    #endif
    for(k = 0; k < Nz-1; k++)
    {
        Hy[t+1 + k*Nt] = Hy[t + k*Nt] - dt/MU_0 * (Ex[t + (k+1)*Nt] - Ex[t + k*Nt])/dz;
    }
#else
    #ifdef OPEN_MP_SPACE
        #pragma omp parallel for private(k)
    #endif
    for(k = 0; k < Nz-1; k++)
    {
        Hy[t+1 + k*Nt] = Hy[t + k*Nt] - dt/MU_0 * (Ex[t + (k+1)*Nt] - Ex[t + k*Nt])/dz;
    }
#endif
}

/**
 * Classical (non-fractional) Ex field update in 1D domain
 * 
 * @param dz spatial step size
 * @param Nz domain size (cells)
 * @param dt time step
 * @param Nt length of simulation (timesteps)
 * @param Ex Ex field array 
 * @param Hy Hy field array
 * @param t current timestep
 * @return 
 */
void ExClassicUpdate(const double dz, const int Nz, const int k_bound, const double dt, const int Nt,
                     double* Ex, const double* Hy, const int t)
{
    int k = 1;

#ifdef MUR_CONDITION
    #ifdef OPEN_MP_SPACE
        #pragma omp parallel for private(k)
    #endif
    for(k = 1; k < Nz-1; k++)
    {
        Ex[t+1 + k*Nt] = Ex[t + k*Nt] - dt/EPS_0* (Hy[t+1 + k*Nt] - Hy[t+1 + (k-1)*Nt]) / dz; // update based on Hy rotation
    }
    // Mur condition - left boundary
    Ex[t+1 + 0*Nt] = Ex[t + 1*Nt] + (C_CONST*dt-dz)/(C_CONST*dt+dz) * (Ex[t+1 + 1*Nt] - Ex[t + 0*Nt]);
    // Mur condition - right boundary
    Ex[t+1 + (Nz-1)*Nt] = Ex[t + (Nz-2)*Nt] + (C_CONST*dt-dz)/(C_CONST*dt+dz) * (Ex[t+1 + (Nz-2)*Nt] - Ex[t + (Nz-1)*Nt]);

#else
    #ifdef OPEN_MP_SPACE
        #pragma omp parallel for private(k)
    #endif
    for(k = 1; k < Nz-1; k++)
    {
        Ex[t+1 + k*Nt] = Ex[t + k*Nt] - dt/EPS_0* (Hy[t+1 + k*Nt] - Hy[t+1 + (k-1)*Nt]) / dz; // update based on Hy rotation
    }
#endif

}

/**
 * simulation
 * 
 * @param dz spatial step size
 * @param Nz domain size (cells)
 * @param dt time step
 * @param Nt length of simulation (timesteps)
 * @param alpha fractional order of derivative
 * @param Ex Ex field array
 * @param Hy Hy field array
 * @return simulation time in seconds
 */
double simulation(const double dz, const int Nz, const double dt, const int Nt,
                  const double alpha,
                  double* Ex, double* Hy,
                  double* Ex_source, const int k_source,
                  int save_result)
{
    static int sim_counter = 0;
    int k_bound = 600;
    // Precalculate GL coefficients
    double* GL_coeff_arr = calloc(Nt, sizeof(double)); // GL[0]=w1, GL[1]=w2, ... 
    GL_coeff_arr[0] = fracGLCoeff(1.0, alpha, 0+1); // GL_1
    for(int n = 1; n < Nt; n++)
        GL_coeff_arr[n] = fracGLCoeff(GL_coeff_arr[n-1], alpha, n+1);

    double start_time, end_time;
    start_time = omp_get_wtime();

    // main time loop
    for (int t = 0; t < Nt-1; t++)
    {

#ifdef FRACTIONAL_SIM
        HyUpdate(dz, Nz, k_bound, dt, Nt, Ex, Hy, t, alpha, GL_coeff_arr);
#else
        HyClassicUpdate(dz, Nz, k_bound, dt, Nt, Ex, Hy, t);
#endif
        Hy[t+1 + (k_source-1)*Nt] = -Hy[t+1 + k_source*Nt]; // Ex field update as if wave travelled in left direction

#ifdef FRACTIONAL_SIM
        ExUpdate(dz, Nz, k_bound, dt, Nt, Ex, Hy, t, alpha, GL_coeff_arr);
#else
        ExClassicUpdate(dz, Nz, k_bound, dt, Nt, Ex, Hy, t);
#endif        
        Ex[t+1 + (k_source)*Nt] += Ex_source[t+1]; // soft source
        Ex[t+1 + (k_source-1)*Nt] = 0.0; // remove left-travelling wave
        Hy[t+1 + (k_source-1)*Nt] = 0.0;
    }

    end_time = omp_get_wtime();
    double sim_time = end_time - start_time;

    if(save_result)
    {
        char filename[128];
        // sprintf(filename, ".\\results\\Ex_%d.bin", sim_counter);
        sprintf(filename, ".\\results\\Ex.bin");
        saveFieldToBinary(filename, Ex, Nz, Nt, dz, dt, alpha);
        sprintf(filename, ".\\results\\Hy.bin");
        saveFieldToBinary(filename, Hy, Nz, Nt, dz, dt, alpha);
    }
    sim_counter += 1;

    return sim_time;
}


/**
 * Save parameters of simulation to text file. Append to existing data
 * 
 * @param filename name of binary file
 * @param dz 
 * @param Lz 
 * @param Nz  
 * @param dt 
 * @param T 
 * @param Nt 
 * @param alpha 
 * @param sim_time simulation time in seconds
 * @return 
 */
void saveSimParamsToTxt(const char *filename,
                           const double dz, const double Lz, const unsigned int Nz,
                           const double dt, const double T, const unsigned int Nt,
                           const double alpha, const double sim_time)
{
    FILE *fptr;
	if ((fptr = fopen(filename, "a")) == NULL)
    {
        printf("cannot open file!\n");
        exit(1); 
    }
    
    // Get simulation flags
    unsigned int sim_flag = 0;
#ifdef FRACTIONAL_SIM
    sim_flag = sim_flag | 1;
#endif
#ifdef MUR_CONDITION
    sim_flag = sim_flag | 1 << 1;
#endif
#ifdef OPEN_MP_SPACE
    sim_flag = sim_flag | 1 << 2;
#endif
#ifdef TIME_ROW_WISE
    sim_flag = sim_flag | 1 << 3;
#endif

    // Write params separated by ,. Every new data pack - new line
    fprintf(fptr, "%d, ", sim_flag);
    fprintf(fptr, "%e, %d, %e, ", dz, Nz, Lz);
    fprintf(fptr, "%e, %d, %e, ", dt, Nt, T);
    fprintf(fptr, "%e, %e\n", alpha, sim_time);

	fclose(fptr);
}


/**
 * Save parameters of simulation to binary file. Append to existing data.
 * Outdated - moved to saving to txt
 * 
 * @param filename name of binary file
 * @param dz 
 * @param Lz 
 * @param Nz  
 * @param dt 
 * @param T 
 * @param Nt 
 * @param alpha 
 * @param sim_time simulation time in seconds
 * @return 
 */
void saveSimParamsToBinary(const char *filename,
                           const double dz, const double Lz, const unsigned int Nz,
                           const double dt, const double T, const unsigned int Nt,
                           const double alpha, const double sim_time)
{
    FILE *fptr;
	
	if ((fptr = fopen(filename, "a")) == NULL)
    {
        printf("cannot open file!\n");
        exit(1); 
    }

    // Write simulation flags
    unsigned int sim_flag = 0;
#ifdef FRACTIONAL_SIM
    sim_flag = sim_flag | 1;
#endif
#ifdef MUR_CONDITION
    sim_flag = sim_flag | 1 << 1;
#endif
#ifdef OPEN_MP_SPACE
    sim_flag = sim_flag | 1 << 2;
#endif
#ifdef TIME_ROW_WISE
    sim_flag = sim_flag | 1 << 3;
#endif

    size_t written_elements = 0;
    written_elements += fwrite(&sim_flag, sizeof(unsigned int), 1, fptr);
    written_elements += fwrite(&dz, sizeof(double), 1, fptr);
    written_elements += fwrite(&Lz, sizeof(double), 1, fptr);
    written_elements += fwrite(&Nz, sizeof(unsigned int), 1, fptr);
    written_elements += fwrite(&dt, sizeof(double), 1, fptr);
    written_elements += fwrite(&T, sizeof(double), 1, fptr);
    written_elements += fwrite(&Nt, sizeof(unsigned int), 1, fptr);
    written_elements += fwrite(&alpha, sizeof(double), 1, fptr);
    written_elements += fwrite(&sim_time, sizeof(double), 1, fptr);

    printf("written_elements: %d\n",written_elements); // TEMP
	fclose(fptr);
}


/** Read params from file. Test to validate saved data.
 * \param filename name of file
 */
void FdtdCpuReadPlanefromFile(const char *filename)
{
	FILE *fptr;
	if ((fptr = fopen(filename, "rb")) == NULL)
	{
		perror("Error");
		printf("FdtdCpuReadPlanefromFile, cannot open file");
        exit(1);
	}
    unsigned int sim_flag;
    double dz;
    double Lz;
    unsigned int Nz;
    double dt;
    double T;
    unsigned int Nt;
    double alpha;
    double sim_time;

    for(int i = 0; i < 50; i++) // hardcoded number of data to read
    {
        printf("i = %d:\n", i);
        fread(&sim_flag, sizeof(unsigned int), 1, fptr);

        fread(&dz, sizeof(double), 1, fptr);
        fread(&Lz, sizeof(unsigned int), 1, fptr);
        fread(&Nz, sizeof(unsigned int), 1, fptr);

        fread(&dt, sizeof(double), 1, fptr);
        fread(&T, sizeof(double), 1, fptr);
        fread(&Nt, sizeof(unsigned int), 1, fptr);

        fread(&alpha, sizeof(double), 1, fptr);
        fread(&sim_time, sizeof(double), 1, fptr);

        printf("sim_flag = %d, dz = %e, Nz = %d, dt = %e, T = %e, Nt = %d, alpha = %e, sim_time = %e\n", 
                sim_flag, dz, Nz, dt, T, Nt, alpha, sim_time);
    }
	fclose(fptr);
}



/**
 * Save field to binary file
 * 
 * @param filename name of binary file
 * @param data pointer to 2D-field array to save
 * @param Nz domain size (cells)
 * @param Nt length of simulation (timesteps)
 * @return  
 */
void saveFieldToBinary(const char *filename,
                        const double *data,
                        const unsigned int Nz,
                        const unsigned int Nt,
                        const double dz,
                        const double dt,
                        const double alpha)
{
    FILE *fptr;
	
	if ((fptr = fopen(filename, "wb")) == NULL)
    {
        printf("cannot open file!\n");
        exit(1); 
    }
		
	// Writing dimensions of simulation domain and step sizes to file
    fwrite(&Nz, sizeof(unsigned int), 1, fptr);
	fwrite(&Nt, sizeof(unsigned int), 1, fptr);
    fwrite(&dz, sizeof(double), 1, fptr);
	fwrite(&dt, sizeof(double), 1, fptr);
    fwrite(&alpha, sizeof(double), 1, fptr);


	fwrite(data, sizeof(double), Nz*Nt, fptr);
	fclose(fptr);
}

/**
 * Calculate next coefficient for Grunwald-Letnikov derivative
 * based on previous coefficient, order of derivative and coefficient number 
 * 
 * @param w previous coefficient
 * @param alpha fractional order of derivative 
 * @param n coefficient number (index)
 * @return new (next) coefficient
 */
double fracGLCoeff(const double w, const double alpha, const int n)
{
    return (1.0 - (1.0+alpha)/n) * w;
}


/**
 * Calculate generalised binomial coefficient
 * 
 * @param alpha fractional order of derivative
 * @param k coefficient number (index)
 * @return binomial coefficient
 */
double binomialCoeff(const double alpha, const int k)
{
    return tgamma(alpha+1)/tgamma(alpha+1-k)/tgamma(k+1);
}


/**
 * Find smallest dt for which simulation is unstable, for different orders (alpha).
 * For every considered order alpha multiple short simulations, with dt values near analytical
 * stability boundary are performed.
 * Simulation is considered unstable if max value of abs(Ex) exceeds certain treshold.
 * 
 */
void checkStability()
{
    double alpha_array[ALPHA_ST_NUM] = {0.995, 0.99, 0.98, 0.97, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.51};
    double unstable_dt_array[ALPHA_ST_NUM] = {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0}; // smallest dt, for which simulation is unstable

    int Nz = 120;
    int Nt =  200;
    int Nt_base = 150;
    int k_source = 60;
    double dz = 0.02e-6;

    // iterate oveer alpha
    for(int i=0; i < ALPHA_ST_NUM; i++)
    {
        double alpha = alpha_array[i];
        double dt_base = pow(2.0, 1.0-1.0/alpha) * pow(sqrt(EPS_0*MU_0) * dz, 1.0/alpha); // analytically calculated boundary
        double dt_array[DT_ST_NUM]; // dt values to check
        // double dt_resolution = 0.02e-17;
        double dt_resolution = dt_base/100.0;

        Nt = (int) (Nt_base*(1/alpha * 1/alpha)); // smaller alpha -> smaller exponential growth in unstable case -> 
                                                  // takes longer to decide weather simulation is unstable

        for(int j = 0; j < DT_ST_NUM; j++)  // [dt-7*dt_r, ..., dt-dt_r]
            dt_array[j] = dt_base - 7*dt_resolution + j*dt_resolution + dt_resolution/2;

        printf("i=%d, Nt=%d\n", i, Nt);

        double unstable_dt = -1.0; // smallest dt, for which simulation is unstable
        // iterate over dt values
        for(int j = 0; j < DT_ST_NUM; j++)
        {
            double dt = dt_array[j];
            double* Ex = calloc(Nz*Nt, sizeof(double));
            double* Hy = calloc(Nz*Nt, sizeof(double));
            double* Ex_source = calloc(Nt, sizeof(double));

            // init condition
            Ex[0+ k_source*Nt] = 1.0;

            // perform short simulation
            simulation(dz, Nz, dt, Nt, alpha, Ex, Hy, Ex_source, k_source, 1);

            // If max_value greater than treshold - there is growth and simulation is unstable
            // Field value can exceed init value for first timesteps, but it decreases later.
            // Thats why greater threshold was chosen
            // double max_value = findMaxAbsValue(Nz, Nt, Ex);
            double max_value = findMaxAbsValueLastTimeStep(Nz, Nt, Ex);
            if(fabs(max_value) > 1.0)
            {
                if(unstable_dt < 0.0)
                    unstable_dt = dt;
            }

            printf("alpha=%f, dt=%e, max_value=%e\n", alpha, dt, max_value);
            free(Ex);
            free(Hy);
            free(Ex_source);
        }
        unstable_dt_array[i] = unstable_dt;
        printf("\n");
    }
    printf("alpha_arr_num= [");
    for(int i=0; i < ALPHA_ST_NUM-1; i++)
        printf("%e, ", alpha_array[i]);
    printf("%e]\n", alpha_array[ALPHA_ST_NUM-1]);

    printf("unstable_dt= [");
    for(int i=0; i < ALPHA_ST_NUM-1; i++)
        printf("%e, ", unstable_dt_array[i]);
    printf("%e]\n", unstable_dt_array[ALPHA_ST_NUM-1]);
}


/**
 * Find maximum absolute value in field array
 * 
 * @param Nz space size of array
 * @param Nt temporal size of array
 * @param Ex field array
 * @return maximum value of filed in given array
 */
double findMaxAbsValue(const int Nz, const int Nt, double* Ex)
{
    double max_value = 0.0;
    int change_counter = 0;
    for(int i=0; i < Nt*Nz; i++)
    {
        if(fabs(Ex[i]) > fabs(max_value))
        {
            max_value = Ex[i];
            change_counter += 1;
        }
    }

    return max_value;
}



/**
 * Find maximum absolute value in field array
 * 
 * @param Nz space size of array
 * @param Nt temporal size of array
 * @param Ex field array
 * @return maximum value of filed in given array
 */
double findMaxAbsValueLastTimeStep(const int Nz, const int Nt, double* Ex)
{
    double max_value = 0.0;
    int change_counter = 0;
    for(int k=0; k < Nz; k++)
    {
        if(fabs(Ex[(Nt-1) + k*Nt]) > fabs(max_value))
        {
            max_value = Ex[(Nt-1) + k*Nt];
            change_counter += 1;
        }
    }

    return max_value;
}