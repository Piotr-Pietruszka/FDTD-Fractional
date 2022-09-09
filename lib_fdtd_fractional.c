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
    #pragma omp parallel for
#endif
    for(k = 0; k < Nz-1; k++)
    {
        // update based on Hy rotation
        Hy[t+1 + k*Nt] =  -pow(dt, alpha)/MU_0 * (Ex[t + (k+1)*Nt] - Ex[t + k*Nt])/dz;

        // update based on previous Hy values (from GL derivative)
        for(n = 0; n < t; n++)
        {
            Hy[t+1 + k*Nt] -= Hy[t-n + k*Nt] * GL_coeff_arr[n];
        }
    }
#else
    int k = 0;
    int n = 0;
#ifdef OPEN_MP_SPACE
    #pragma omp parallel for
#endif
    for(k = 0; k < Nz-1; k++)
    {
        // update based on Ex rotation
        Hy[(t+1)*Nz + k] =  -pow(dt, alpha)/MU_0 * (Ex[t*Nz + k+1] - Ex[t*Nz + k])/dz;

        // update based on previous Hy values (from GL derivative)
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
#ifdef OPEN_MP_SPACE
    #pragma omp parallel for
#endif
    for(k = 1; k < Nz-1; k++)
    {
        // update based on Hy rotation
        Ex[t+1 + k*Nt] = -pow(dt, alpha)/EPS_0* (Hy[t+1 + k*Nt] - Hy[t+1 + (k-1)*Nt]) / dz; // update based on Hy rotation
        
        // update based on previous Ex values (from GL derivative)
        for(n = 0; n < t; n++)
        {
            Ex[t+1 + k*Nt] -= Ex[t-n + k*Nt] * GL_coeff_arr[n];
        }
    }
#ifdef MUR_CONDITION  
    // Mur condition - left boundary
    Ex[t+1 + 0*Nt] = (C_CONST*pow(dt, alpha)-dz)/(C_CONST*pow(dt, alpha)+dz) * Ex[t+1 + 1*Nt] +
                     (C_CONST*pow(dt, alpha))/(C_CONST*pow(dt, alpha)+dz) * (Ex[t + 1*Nt] - Ex[t + 0*Nt]);
    for(n = 0; n < t; n++)
    {
        Ex[t+1 + 0*Nt] -=  (dz)/(C_CONST*pow(dt, alpha)+dz) * GL_coeff_arr[n] * Ex[t-n + 1*Nt];
        Ex[t+1 + 0*Nt] -=  (dz)/(C_CONST*pow(dt, alpha)+dz) * GL_coeff_arr[n] * Ex[t-n + 0*Nt];
    }
    // Mur condition - right boundary   
    Ex[t+1 + (Nz-1)*Nt] = (C_CONST*pow(dt, alpha)-dz)/(C_CONST*pow(dt, alpha)+dz) * Ex[t+1 + (Nz-2)*Nt] +
                          (C_CONST*pow(dt, alpha))/(C_CONST*pow(dt, alpha)+dz) * (Ex[t + (Nz-2)*Nt] - Ex[t + (Nz-1)*Nt]);
    for(n = 0; n < t; n++)
    {
        Ex[t+1 + (Nz-1)*Nt] -=  (dz)/(C_CONST*pow(dt, alpha)+dz) * GL_coeff_arr[n] * Ex[t-n + (Nz-1)*Nt];
        Ex[t+1 + (Nz-1)*Nt] -=  (dz)/(C_CONST*pow(dt, alpha)+dz) * GL_coeff_arr[n] * Ex[t-n + (Nz-2)*Nt];
    }
#endif
#else
    int k = 1;
    int n = 0;
#ifdef OPEN_MP_SPACE
    #pragma omp parallel for
#endif
    for(k = 1; k < Nz-1; k++)
    {
        // update based on Hy rotation
        Ex[(t+1)*Nz + k] = -pow(dt, alpha)/EPS_0* (Hy[(t+1)*Nz + k] - Hy[(t+1)*Nz + k-1]) / dz; // update based on Hy rotation
        
        // update based on previous Ex values (from GL derivative)
        for(n = 0; n < t; n++)
        {
            Ex[(t+1)*Nz + k] -= Ex[(t-n)*Nz + k] * GL_coeff_arr[n];
        }
    }
#ifdef MUR_CONDITION
    // Mur condition - left boundary
    Ex[t+1 + 0*Nt] = (C_CONST*pow(dt, alpha)-dz)/(C_CONST*pow(dt, alpha)+dz) * Ex[(t+1)*Nz + 1] +
                     (C_CONST*pow(dt, alpha))/(C_CONST*pow(dt, alpha)+dz) * (Ex[(t)*Nz + 1] - Ex[(t)*Nz + 0*Nt]);
    for(n = 0; n < t; n++)
    {
        Ex[(t+1)*Nz + 0] -=  (dz)/(C_CONST*pow(dt, alpha)+dz) * GL_coeff_arr[n] * Ex[(t-n)*Nz + 1];
        Ex[(t+1)*Nz + 0] -=  (dz)/(C_CONST*pow(dt, alpha)+dz) * GL_coeff_arr[n] * Ex[(t-n)*Nz + 0];
    }
    // Mur condition - right boundary   
    Ex[(t+1)*Nz + Nz-1] = (C_CONST*pow(dt, alpha)-dz)/(C_CONST*pow(dt, alpha)+dz) * Ex[(t+1)*Nz + Nz-2] +
                          (C_CONST*pow(dt, alpha))/(C_CONST*pow(dt, alpha)+dz) * (Ex[(t)*Nz + Nz-2] - Ex[(t)*Nz + Nz-1]);
    for(n = 0; n < t; n++)
    {
        Ex[(t+1)*Nz + Nz-1] -=  (dz)/(C_CONST*pow(dt, alpha)+dz) * GL_coeff_arr[n] * Ex[(t-n)*Nz + Nz-1];
        Ex[(t+1)*Nz + Nz-1] -=  (dz)/(C_CONST*pow(dt, alpha)+dz) * GL_coeff_arr[n] * Ex[(t-n)*Nz + Nz-2];
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
#ifdef OPEN_MP_SPACE
    #pragma omp parallel for
#endif
    for(k = 0; k < Nz-1; k++)
    {
        Hy[t+1 + k*Nt] = Hy[t + k*Nt] - dt/MU_0 * (Ex[t + (k+1)*Nt] - Ex[t + k*Nt])/dz;
    }
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

#ifdef OPEN_MP_SPACE
    #pragma omp parallel for
#endif
    for(k = 1; k < Nz-1; k++)
    {
        Ex[t+1 + k*Nt] = Ex[t + k*Nt] - dt/EPS_0* (Hy[t+1 + k*Nt] - Hy[t+1 + (k-1)*Nt]) / dz; // update based on Hy rotation
    }
#ifdef MUR_CONDITION
    // Mur condition - left boundary
    Ex[t+1 + 0*Nt] = Ex[t + 1*Nt] + (C_CONST*dt-dz)/(C_CONST*dt+dz) * (Ex[t+1 + 1*Nt] - Ex[t + 0*Nt]);
    // Mur condition - right boundary
    Ex[t+1 + (Nz-1)*Nt] = Ex[t + (Nz-2)*Nt] + (C_CONST*dt-dz)/(C_CONST*dt+dz) * (Ex[t+1 + (Nz-2)*Nt] - Ex[t + (Nz-1)*Nt]);
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

#ifdef ADD_SOURCE
        Hy[t+1 + (k_source-1)*Nt] = -Hy[t+1 + k_source*Nt]; // Ex field update as if wave travelled in left direction
#endif

#ifdef FRACTIONAL_SIM
        ExUpdate(dz, Nz, k_bound, dt, Nt, Ex, Hy, t, alpha, GL_coeff_arr);
#else
        ExClassicUpdate(dz, Nz, k_bound, dt, Nt, Ex, Hy, t);
#endif

#ifdef ADD_SOURCE
        Ex[t+1 + (k_source)*Nt] += Ex_source[t+1]; // soft source
        Ex[t+1 + (k_source-1)*Nt] = 0.0; // remove left-travelling wave
        Hy[t+1 + (k_source-1)*Nt] = 0.0;
#endif
    }

    end_time = omp_get_wtime();
    double sim_time = end_time - start_time;
    printf("Simulation finished. Saving results to files\n");
    if(save_result)
    {
        char filename[128];
        sprintf(filename, ".\\results\\Ex.bin");
        saveFieldToBinary(filename, Ex, Nz, Nt, dz, dt, alpha);
        printf("Ex field saved\n");
        sprintf(filename, ".\\results\\Hy.bin");
        saveFieldToBinary(filename, Hy, Nz, Nt, dz, dt, alpha);
        printf("Hy field saved\n");
    }
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

    // Writing data in chunks - for data larger than 4 GB at once fwrite can hang
    unsigned int offset = 0;
    unsigned int chunk_size = CHUNK_SIZE_BYTES / sizeof(double); // chunk size - 2 GB equivalent
    if(Nz*Nt > chunk_size)
    {
        // printf("Writing in chunks. Single chunk size: %d\n", chunk_size);
        unsigned int last_data = 0;
        for (unsigned int i = 0; i < Nz*Nt-chunk_size; i += chunk_size)
        {
            int ret = fwrite(data+i, sizeof(double), chunk_size, fptr);
            if (ret != chunk_size) {
                printf("Stream error indication %d\n", ferror(fptr));
            }
            last_data = i+chunk_size;
        }
        int ret = fwrite(data+last_data, sizeof(double), Nz*Nt-last_data, fptr);
        if (ret != Nz*Nt-last_data) {
            printf("Stream error indication %d\n", ferror(fptr));
        }
        
    }
    else
    {
        int ret = fwrite(data, sizeof(double), Nz*Nt, fptr);
        if (ret != Nz*Nt) {
            printf("Stream error indication %d\n", ferror(fptr));
        }
    }
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
