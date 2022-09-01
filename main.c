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
#else
    double alpha = 1.0;
    double dt = 0.999*dz/C_CONST; 
#endif
    enum SourceType source_type = SINUS;
    // enum SourceType source_type = TRIANGLE;


    // double Lz = 140.0e-6;
    double Lz = 80.0e-6;
    double T = 18e-14;
    // double T = 20e-14;
    char source_char = 'm';
/*
    printf("Simulation of electromagnetic wave propagation in fractional-order material\n");
    int correct_input = 0;
    while(!correct_input)
    {
        // Get alpha from user
        printf("\nOrder - alpha [0.9 - 1.0]: ");
        if (!scanf("%lf", &alpha))
        {
            scanf("%*[^\n]"); //discard that line up to the newline
            printf("Wrong value of order!\n");
            continue;
        }
        else if(alpha < 0.9 || alpha > 1.0)
        {
            printf("Wrong value of order - use values only from [0.9, 1] interval\n");
            continue;
        }

        // Get Lz from user
        printf("Size of computational domain - Lz [1e-6 - 200e-6]: ");
        if (!scanf("%lf", &Lz))
        {
            scanf("%*[^\n]"); //discard that line up to the newline
            printf("Wrong value of Lz\n");
            continue;
        }
        else if(Lz < 1e-6 || Lz > 200e-6)
        {
            printf("Wrong value of Lz - use values only from [1e-6, 200e-6] interval\n");
            continue;
        }
        
        // Get T from user
        printf("Simulation time - T [1e-14 - 40e-14]: ");
        if (!scanf("%lf", &T))
        {
            scanf("%*[^\n]"); //discard that line up to the newline
            printf("Wrong value of T!\n");
            continue;
        }
        else if(T < 1e-14 || T > 40e-14)
        {
            printf("Wrong value of T - use values only from [1e-14 - 40e-14] interval\n");
            continue;
        }  
        scanf("%*[^\n]"); //discard that line up to the newline
        // Get source type
        printf("Source type: m (modulated gaussian) / r (rectangle) / t (triangle) / g (gaussian): ");
        if (!scanf(" %c", &source_char))
        {
            scanf("%*[^\n]"); //discard that line up to the newline
            printf("Wrong value of source type!\n");
            continue;
        }
        else if(!(source_char=='m' || source_char=='r' || source_char=='t' || source_char=='g'))
        {
            printf("Wrong value of source_type - use values only m, r, t or\n");
            continue;
        }  

        correct_input = 1;
    }
    */
    
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

    if(Ex==NULL || Hy == NULL || Ex_source == NULL) 
    {
        printf("Error - couldn't allocate memory for arrays\n"); 
        exit(1); 
    } 

    // Source
    int k_source = (int) (0.1e-6/dz);
    printf("k_source= %d\n", k_source);


    // First and last timesteps of source signal
    int first_t_source = -1;
    int last_t_source = -1;
    int stop_sin = -1;

    for (int t = 0; t < Nt-1; t++)
    {
        
        if (source_type == MODULATED_GAUSSIAN)
        {
            // double delay = 7.957747154594768e-15 - dt/2.0;
            // Ex_source[t] = 5.9916e+08 * pow(dt, alpha) / dz * cos((t*dt-delay)*3.707079331235956e+15) * exp( -pow((t*dt-delay) / (1.989436788648692e-15), 2.0) ); // modulated gaussian       
            
            double max_fr = 750e12;
            double min_fr = 430e12;
            double central_fr = (max_fr+min_fr) / 2.0;
            double w_central = 2* M_PI * central_fr;

            double tau = 2.0 / M_PI / (max_fr-min_fr);
            double delay = 4*tau+dt/2;

            Ex_source[t] = 5.9916e+08 * pow(dt, alpha) / dz * cos((t*dt-delay)*w_central) * exp( -pow((t*dt-delay) / tau, 2.0) ); // modulated gaussian

            // // double delay = 7.957747154594768e-15 - dt/2.0;
            // Ex_source[t] = 5.9916e+08 * pow(dt, alpha) / dz * cos((t*dt-delay)*3.707079331235956e+15) * exp( -pow((t*dt-delay) / (1.989436788648692e-15), 2.0) ); // modulated gaussian
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
        else if(source_type == SINUS)
        {
            double w_sin = 2* M_PI * 3e14;
            if(t*dt > 8e-14 && sin((t*dt)*w_sin) < 1e-4)
                stop_sin = 1;
            if(stop_sin < 0)
            {
                Ex_source[t] = sin((t*dt)*w_sin); // modulated gaussian       
            }        
        }
    }
    if(source_type == RECTANGLE)
    {
        Ex_source[first_t_source-1] = 1.5; Ex_source[first_t_source-2] = 1.0; Ex_source[first_t_source-3] = 0.5; 
        Ex_source[last_t_source+1] = 1.5; Ex_source[last_t_source+2] = 1.0; Ex_source[last_t_source+3] = 0.5;     
    
        // Ex_source[first_t_source-1] = 1.75; Ex_source[first_t_source-2] = 1.5; Ex_source[first_t_source-3] = 1.25; Ex_source[first_t_source-4] = 1; Ex_source[first_t_source-5] = 0.75; Ex_source[first_t_source-6] = 0.5; Ex_source[first_t_source-7] = 0.25;
        // Ex_source[last_t_source+1] = 1.75; Ex_source[last_t_source+2] = 1.5; Ex_source[last_t_source+3] = 1.25; Ex_source[last_t_source+4] = 1; Ex_source[last_t_source+5] = 0.75; Ex_source[last_t_source+6] = 0.5; Ex_source[last_t_source+7] = 0.25;    
    }

    char filename[128];
    sprintf(filename, ".\\results\\source.bin");
    saveFieldToBinary(filename, Ex_source, 1, Nt, dz, dt, alpha);

    double sim_time = simulation(dz, Nz, dt, Nt, alpha, Ex, Hy, Ex_source, k_source, 1);
    printf("simulation time: %lf s\n", sim_time);

    sprintf(filename, ".\\results\\results.txt");
    saveSimParamsToTxt(filename, dz, Lz, Nz,
                        dt, T, Nt, alpha, sim_time);

    // free allocated memory
    free(Ex);
    free(Hy);
    free(Ex_source);
}