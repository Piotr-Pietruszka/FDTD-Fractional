# Fractional simulation
Simulation of electromagnetic wave propagation in fractional order material, using FDTD method. Simulation is written in C and visualisation of results in MATLAB.

##  Execution:
Compile program, by running in project folder: 

    gcc -g *.c -o main.exe -O3 -Wall -fopenmp

Run: 

    main.exe


## Files
### Simulation:
- lib_fdtd_fractional.h - header file with function declarations, defined simulation type and physical constants.
- lib_fdtd_fractional.c - function's definitions: update of field values, source initialization, saving results to binary file, simulation time loop.
- main.c - main file with defined simulation parameters and function calls.

### Visualisation - Matlab scripts:
- toWorkspace.m - load simulation parameters to workspace (order of fractional material, space and time step size, etc).
- movie.m - plot field values in space domain, for different time steps. Save result as .avi file. Example use: `movie("Ex.bin");`. 
- plotFieldAtPositions.m - plot field value in time, at given positions. Example use: `plotFieldAtPositions('Ex.bin', [10, 4010]);`. 
