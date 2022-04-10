#include "lib_fdtd_fractional.h"

/**
 * Update Hy field in 1D domain
 * 
 * @param dz - spatial step size
 * @param Nz domain size (cells)
 * @param dt time step
 * @param Nt length of simulation (timesteps)
 * @param Ex Ex field array 
 * @param Hy Hy field array
 * @param t current timestep
 * @param alpha fractional order of derivative
 * @return 
 */
void HyUpdate(const double dz, const int Nz, const double dt, const int Nt,
                const double* Ex, double* Hy, const int t, const double alpha)
{
    int k = 0;
    int n = 0;
    
#pragma omp parallel for private(k)
    for(k = 0; k < Nz-1; k++)
    {
        // update based on Hy rotation
        Hy[(t+1)*Nz + k] =  -pow(dt, alpha)/MU_0 * (Ex[t*Nz + k+1] - Ex[t*Nz + k])/dz;

        // update based on previous Hy values (from GL derivative)
        double w = 1.0;
        for(n = 0; n < t; n++)
        {
            w = fracGLCoeff(w, alpha, n+1);
            Hy[(t+1)*Nz + k] -= Hy[(t-n)*Nz + k] * w;
            // Hy[(t+1)*Nz + k] += Hy[(t-n)*Nz + k] * pow(-1.0, n+1) * binomialCoeff(alpha, n+1);
        }
    }
}


/**
 * Update Ex field in 1D domain
 * 
 * @param dz - spatial step size
 * @param Nz domain size (cells)
 * @param dt time step
 * @param Nt length of simulation (timesteps)
 * @param Ex Ex field array 
 * @param Hy Hy field array
 * @param t current timestep
 * @param alpha fractional order of derivative
 * @return 
 */
void ExUpdate(const double dz, const int Nz, const double dt, const int Nt,
                double* Ex, const double* Hy, const int t, const double alpha)
{
    int k = 1;
    int n = 0;
#pragma omp parallel for private(k)
    for(k = 1; k < Nz; k++)
    {
        // update based on Hy rotation
        Ex[(t+1)*Nz + k] = -pow(dt, alpha)/EPS_0* (Hy[(t+1)*Nz + k] - Hy[(t+1)*Nz + k-1]) / dz; // update based on Hy rotation
        
        // update based on previous Ex values (from GL derivative)
        double w = 1.0;
        for(n = 0; n < t; n++)
        {
            w = fracGLCoeff(w, alpha, n+1);
            // Ex[(t+1)*Nz + k] += Ex[(t-n)*Nz + k] * pow(-1.0, n+1) * binomialCoeff(alpha, n+1);
            Ex[(t+1)*Nz + k] -= Ex[(t-n)*Nz + k] * w;
        }      
    }
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
 * @return 
 */
void simulation(const double dz, const int Nz, const double dt, const int Nt,
                const double alpha,
                double* Ex, double* Hy,
                double* Ex_inc, double* Hy_inc)
{
    int k_source = 300;
    // main time loop
    for (int t = 0; t < Nt-1; t++)
    {
        Ex_inc[t] = sin(t*dt*2*3.14/0.15e-14) * exp( -pow((t*dt-0.75e-14) / (0.2e-14), 2.0) );

        HyUpdate(dz, Nz, dt, Nt, Ex, Hy, t, alpha);
        // TODO: tfsf Hy update
        // Hy[(t+1)*Nz + k_source - 1] += pow(dt, alpha)/MU_0/dz * Ex_inc[t];
        // Hy[(t+1)*Nz + k_source - 1] = 0.0;

        ExUpdate(dz, Nz, dt, Nt, Ex, Hy, t, alpha);
        // TODO: tfsf Ex update
        // Ex[(t+1)*Nz + k_source - 1] += pow(dt, alpha)/EPS_0/dz * Hy_inc[t];
        // Ex[(t+1)*Nz + k_source] = Ex_inc[t];

        // double update_value = sin(t*dt*2*3.14/0.15e-14) * exp( -pow((t*dt-0.75e-14) / (0.2e-14), 2.0) );
        double update_value = sin(t*dt*2*3.14/0.3e-14) * exp( -pow((t*dt-0.75e-14) / (0.2e-14), 2.0) );
        Ex[t*Nz + k_source] += update_value; // soft source
    }

    char filename[128];
    sprintf(filename, ".\\results\\Ex.bin");
    saveFieldToBinary(filename, Ex, Nz, Nt);
}


/**
 * Save filed to binary file
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
                        const unsigned int Nt)
{
    FILE *fptr;
	
	if ((fptr = fopen(filename, "wb")) == NULL)
    {
        printf("cannot open file!\n");
        exit(1); 
    }
		
	// Writing dimensions of simulation domain to file
	fwrite(&Nz, sizeof(unsigned int), 1, fptr);
	fwrite(&Nt, sizeof(unsigned int), 1, fptr);

	fwrite(data, sizeof(double), Nz*Nt, fptr);
	fclose(fptr);
}


/**
 * Calculate next coefficient for Grunwald-Letnikov derivative
 * based on previous coefficient, order of derivative and  
 * 
 * @param w previous coefficient
 * @param alpha fractional order of derivative 
 * @param n 
 * @return new coefficient
 */
double fracGLCoeff(const double w, const double alpha, const int n)
{
    return (1.0 - (1.0+alpha)/n) * w;
}


/**
 * Calculate generalised binomial coefficient
 * 
 * @param alpha fractinal order of derivative
 * @param k
 * @return binomial coefficient
 */
double binomialCoeff(const double alpha, const int k)
{
    return tgamma(alpha+1)/tgamma(alpha+1-k)/tgamma(k+1);
}
