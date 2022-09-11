#include "lib_fdtd_fractional.h"

/**
 * Soft source - function to initialize source field in time (Ex_source). 
 * 
 * @param dz spatial step size
 * @param dt time step size
 * @param Nt length of simulation (timesteps)
 * @param alpha fractional order of derivative
 * @param Ex_source source Ex field array
 * @param source_type type of source signal
 * @return 
 */
void initializeSource(const double dz, const double dt, const int Nt, const double alpha,
                      double* Ex_source, enum SourceType source_type)
{
    int first_t_source = -1;
    int last_t_source = -1;
    for (int t = 0; t < Nt-1; t++)
    {
        if (source_type == MODULATED_GAUSSIAN)
        {
            double max_fr = 750e12;
            double min_fr = 430e12;
            double central_fr = (max_fr+min_fr) / 2.0;
            double w_central = 2* M_PI * central_fr;
            double tau = 2.0 / M_PI / (max_fr-min_fr);
            double delay = 4*tau+dt/2;
            Ex_source[t] = 5.9916e+08 * pow(dt, alpha) / dz * cos((t*dt-delay)*w_central) * 
                           exp( -pow((t*dt-delay) / tau, 2.0) ); // modulated gaussian
        }
        else if (source_type == TRIANGLE)
        {
            if(t*dt > 0.4e-14 && t*dt < 0.5e-14)
                Ex_source[t] = (t*dt - 0.4e-14) * 1 / 0.1e-14; // linear function growing to reach 1 at 0.5
            else if(t*dt > 0.5e-14 && t*dt < 0.6e-14)
                Ex_source[t] = 1.0 - (t*dt - 0.5e-14) * 1 / 0.1e-14; // linear function decreasing to reach 0 at 0.6
        }
        else if(source_type == RECTANGLE)
        {
            double start = 3.387e-15;
            double end = 6.8e-15;
            if(t*dt > start && t*dt < end)
            {
                if(first_t_source < 0)
                    first_t_source = t;
                Ex_source[t] = 2.0; // rectangular function
                last_t_source = t;
            }
        }
        else if(source_type == GAUSSIAN)
        {
            double delay = 7.957747154594768e-15 - dt/2.0;
            Ex_source[t] = 5.9916e+08 * pow(dt, alpha) / dz * exp( -pow((t*dt-delay) / (1.989436788648692e-15), 2.0) ); // modulated gaussian       
        }
    }
    if(source_type == RECTANGLE)
    {
        Ex_source[first_t_source-1] = 1.5; Ex_source[first_t_source-2] = 1.0; Ex_source[first_t_source-3] = 0.5; 
        Ex_source[last_t_source+1] = 1.5; Ex_source[last_t_source+2] = 1.0; Ex_source[last_t_source+3] = 0.5;     
    }
}


/**
 * Total field / scattered field. Function to initialize Ex incident field (Ex_source)
 * and Hy incident field (Hy_source). Incident Hy field is for position (k_b-1/2)*dz and 
 * incident Ex field for k_b*dz
 * 
 * @param dt time step size
 * @param Nt length of simulation (timesteps)
 * @param dt_cl_ideal ideal time-step for vacuum (dz/c)
 * @param Ex_source Ex incident field
 * @param Hy_source EHyx incident field
 * @param source_type type of source signal
 * @return 
 */
void initializeSourceTfsf(const double dt, const int Nt, const double dt_cl_ideal,
                          double* Ex_source,  double* Hy_source, enum SourceType source_type)
{
    if(source_type != MODULATED_GAUSSIAN &&  source_type != GAUSSIAN)
    {
        printf("For TFSF interface: only modulated gaussian and gaussian source waves"
               " are available. Using modulated gaussian.\n");
        source_type = MODULATED_GAUSSIAN;      
    }
    for (int t = 0; t < Nt-1; t++)
    {
        if (source_type == MODULATED_GAUSSIAN)
        {
            double max_fr = 750e12;
            double min_fr = 430e12;
            double central_fr = (max_fr+min_fr) / 2.0;
            double w_central = 2* M_PI * central_fr;
            double tau = 2.0 / M_PI / (max_fr-min_fr);
            double delay = 4*tau+dt/2;
            Hy_source[t] = -sqrt(EPS_0/MU_0) * cos(((t+0.5+0.5*dt_cl_ideal/dt)*dt-delay)*w_central) * exp( -pow(((t+0.5+0.5*dt_cl_ideal/dt)*dt-delay) / tau, 2.0) );
            Ex_source[t] = 1.0 * cos((t*dt-delay)*w_central) * exp( -pow((t*dt-delay) / tau, 2.0) );
        }
        else if(source_type == GAUSSIAN)
        {
            double delay = 7.957747154594768e-15 - dt/2.0;
            Hy_source[t] = -sqrt(EPS_0/MU_0) * exp( -pow(((t+0.5+0.5*dt_cl_ideal/dt)*dt-delay) / (1.989436788648692e-15), 2.0) ); // modulated gaussian
            Ex_source[t] = exp( -pow((t*dt-delay) / (1.989436788648692e-15), 2.0) ); // modulated gaussian       
        }
    }

}


/**
 * Update Hy field in 1D for domain filled with 
 * fractional order material 
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
void HyUpdate(const double dz, const int Nz, const double dt, const int Nt,
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
 #else // TIME_ROW_WISE
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
 #endif // TIME_ROW_WISE
}


/**
 * Update Ex field in 1D domain for domain filled with 
 * fractional order material
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
void ExUpdate(const double dz, const int Nz, const double dt, const int Nt,
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
        Ex[t+1 + k*Nt] = -pow(dt, alpha)/EPS_0* (Hy[t+1 + k*Nt] - Hy[t+1 + (k-1)*Nt]) / dz;

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
  #endif // MUR_CONDITION
 #else // TIME_ROW_WISE
    int k = 1;
    int n = 0;
  #ifdef OPEN_MP_SPACE
    #pragma omp parallel for
  #endif
    for(k = 1; k < Nz-1; k++)
    {
        // update based on Hy rotation
        Ex[(t+1)*Nz + k] = -pow(dt, alpha)/EPS_0* (Hy[(t+1)*Nz + k] - Hy[(t+1)*Nz + k-1]) / dz;

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
  #endif // MUR_CONDITION
 #endif // TIME_ROW_WISE
}


/**
 * Update Hy field in 1D for domain consisting of
 * vaccum (left side) and fractional order material (right side)
 * 
 * @param dz spatial step size
 * @param Nz domain size (cells)
 * @param k_bound index of boundary between materials (vacuum and fr order)
 * @param dt time step
 * @param Nt length of simulation (timesteps)
 * @param Ex Ex field array 
 * @param Hy Hy field array
 * @param t current timestep
 * @param alpha fractional order of derivative
 * @param GL_coeff_arr Nt-size array of Gr-Let derivative coefficients, starts at w1
 * @return 
 */
void HyUpdateDiffMaterials(const double dz, const int Nz, const int k_bound, 
                           const double dt, const int Nt,
                           const double* Ex, double* Hy, const int t, const double alpha,
                           const double* GL_coeff_arr)
{
    int k = 0;
    int n = 0;
 #ifdef OPEN_MP_SPACE
    #pragma omp parallel for
 #endif
    for(k = 0; k <= k_bound; k++) // Vacuum
    {
        Hy[t+1 + k*Nt] = Hy[t + k*Nt] - dt/MU_0 * (Ex[t + (k+1)*Nt] - Ex[t + k*Nt])/dz;
    }
 #ifdef OPEN_MP_SPACE
    #pragma omp parallel for
 #endif
    for(k = k_bound+1; k < Nz-1; k++) // Fractional material
    {
        // update based on Hy rotation
        Hy[t+1 + k*Nt] =  -pow(dt, alpha)/MU_0 * (Ex[t + (k+1)*Nt] - Ex[t + k*Nt])/dz;

        // update based on previous Hy values (from GL derivative)
        for(n = 0; n < t; n++)
        {
            Hy[t+1 + k*Nt] -= Hy[t-n + k*Nt] * GL_coeff_arr[n];
        }
    }
}


/**
 * Update Ex field in 1D domain for domain consisting of
 * vaccum (left side) and fractional order material (right side)
 * 
 * @param dz spatial step size
 * @param Nz domain size (cells)
 * @param k_bound index of boundary between vaccum and fr material 
 * @param dt time step
 * @param Nt length of simulation (timesteps)
 * @param Ex Ex field array 
 * @param Hy Hy field array
 * @param t current timestep
 * @param alpha fractional order of derivative
 * @param GL_coeff_arr Nt-size array of Gr-Let derivative coefficients, starts at w1
 * @return 
 */
void ExUpdateDifferentMaterials(const double dz, const int Nz, const int k_bound, 
                                const double dt, const int Nt,
                                double* Ex, const double* Hy, const int t, const double alpha,
                                const double* GL_coeff_arr)
{
    int k = 1;
    int n = 0;
 #ifdef OPEN_MP_SPACE
    #pragma omp parallel for
 #endif
    for(k = 1; k <= k_bound; k++) // Vacuum
    {
        Ex[t+1 + k*Nt] = Ex[t + k*Nt] - dt/EPS_0* (Hy[t+1 + k*Nt] - Hy[t+1 + (k-1)*Nt]) / dz; // update based on Hy rotation
    }
 #ifdef OPEN_MP_SPACE
    #pragma omp parallel for
 #endif
    for(k = k_bound+1; k < Nz-1; k++)
    {
        // update based on Hy rotation
        Ex[t+1 + k*Nt] = -pow(dt, alpha)/EPS_0* (Hy[t+1 + k*Nt] - Hy[t+1 + (k-1)*Nt]) / dz;

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
 #endif // MUR_CONDITION
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
void HyUpdateClassical(const double dz, const int Nz, const double dt, const int Nt,
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
void ExUpdateClassical(const double dz, const int Nz, const double dt, const int Nt,
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
                  double* Ex_source, double* Hy_source, const int k_source, 
                  const int k_bound, int save_result)
{
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
        // Update Hy
 #if SIMULATION_TYPE == FRACTIONAL_SIMULATION
        HyUpdate(dz, Nz, dt, Nt, Ex, Hy, t, alpha, GL_coeff_arr);
 #elif SIMULATION_TYPE == DIFFERENT_MATERIALS
        HyUpdateDiffMaterials(dz, Nz, k_bound, dt, Nt, Ex, Hy, t, alpha, GL_coeff_arr);
 #elif SIMULATION_TYPE == CLASSICAL_SIMULATION
        HyUpdateClassical(dz, Nz, dt, Nt, Ex, Hy, t);
 #endif
        // Source - Hy modification
 #if SIMULATION_TYPE == DIFFERENT_MATERIALS
        Hy[t+1 + (k_source-1)*Nt] += dt/(MU_0*dz)*Ex_source[t+1]; // TFSF
 #else
        Hy[t+1 + (k_source-1)*Nt] = -Hy[t+1 + k_source*Nt]; // Ex field update as if wave travelled in left direction
 #endif

        // Update Ex
 #if SIMULATION_TYPE == FRACTIONAL_SIMULATION
        ExUpdate(dz, Nz, dt, Nt, Ex, Hy, t, alpha, GL_coeff_arr);
 #elif SIMULATION_TYPE == DIFFERENT_MATERIALS
        ExUpdateDifferentMaterials(dz, Nz, k_bound, dt, Nt, Ex, Hy, t, alpha, GL_coeff_arr);
 #elif SIMULATION_TYPE == CLASSICAL_SIMULATION
        ExUpdateClassical(dz, Nz, dt, Nt, Ex, Hy, t);
 #endif

        // Source - Ex modification
 #if SIMULATION_TYPE == DIFFERENT_MATERIALS
        Ex[t+1 + (k_source)*Nt] -= dt/(EPS_0*dz)*Hy_source[t+1]; // TFSF
 #else
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
        saveFieldToBinary(filename, Ex, Nz, Nt, dz, dt, alpha, k_bound);
        printf("Ex field saved\n");
        sprintf(filename, ".\\results\\Hy.bin");
        saveFieldToBinary(filename, Hy, Nz, Nt, dz, dt, alpha, k_bound);
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
#if SIM_TYPE == FRACTIONAL_SIMULATION
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
                        const double alpha, 
                        const int k_bound)
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
    fwrite(&k_bound, sizeof(int), 1, fptr);

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
