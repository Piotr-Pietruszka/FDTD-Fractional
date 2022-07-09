
filename = "results.bin";
fptr=fopen(filename);

FRACTIONAL_SIM_arr = [];
MUR_CONDITION_arr = [];
OPEN_MP_SPACE_arr = [];
TIME_ROW_WISE_arr = [];
dz_arr = [];
Lz_arr = [];
Nz_arr = [];
dt_arr = [];
T_arr = [];
Nt_arr = [];
alpha_arr = [];
sim_time_arr = [];
sim_flags_arr = [];

%% Read siumlation paramaters for all simulations from file
i = 0;
sim_flags = fread(fptr, 1, 'uint');
while sim_flags
    i = i + 1;
    
    %% Read paramteres for single simulation from file
    FRACTIONAL_SIM = boolean(bitand(sim_flags, 1));
    MUR_CONDITION = boolean(bitand(sim_flags, 2));
    OPEN_MP_SPACE = boolean(bitand(sim_flags, 4));
    TIME_ROW_WISE = boolean(bitand(sim_flags, 8));

    dz = fread(fptr, 1, 'double');
    Lz = fread(fptr, 1, 'double');
    Nz = fread(fptr, 1, 'uint');

    dt = fread(fptr, 1, 'double');
    T = fread(fptr, 1, 'double');
    Nt = fread(fptr, 1, 'uint');

    alpha = fread(fptr, 1, 'double');
    sim_time = fread(fptr, 1, 'double');
    
    %% save to arrays 
    sim_flags_arr(i) = sim_flags;
    FRACTIONAL_SIM_arr(i) = FRACTIONAL_SIM; MUR_CONDITION_arr(i) = MUR_CONDITION; 
    OPEN_MP_SPACE_arr(i) = OPEN_MP_SPACE; TIME_ROW_WISE_arr(i) = TIME_ROW_WISE;
    dz_arr(i) = dz; Lz_arr(i) = Lz; Nz_arr(i) = Nz;
    dt_arr(i) = dt; T_arr(i) = T; Nt_arr(i) = Nt;
    alpha_arr(i) = alpha;
    sim_time_arr(i) = sim_time;


    %% for next iteration
    sim_flags = fread(fptr, 1, 'uint');
    
end
fclose(fptr);

%%
