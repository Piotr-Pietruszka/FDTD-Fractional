
%% Read data from file
filename = "results.txt";
d=importdata("results.txt");




%% Values from data
% Flags
sim_flags_arr = d(:, 1);
FRACTIONAL_SIM_arr = boolean(bitand(sim_flags_arr, 1));
MUR_CONDITION_arr = boolean(bitand(sim_flags_arr, 2));
OPEN_MP_SPACE_arr = boolean(bitand(sim_flags_arr, 4));
TIME_ROW_WISE_arr = boolean(bitand(sim_flags_arr, 8));

% space dims
dz_arr = d(:, 2);
Nz_arr = d(:, 3);
Lz_arr = d(:, 4);

% time dims
dt_arr = d(:, 5);
Nt_arr = d(:, 6);
T_arr = d(:, 7);

% Other - order alpha adn time of sim
alpha_arr = d(:, 8);
sim_time_arr = d(:, 9);


%% Ids to filter by - flags

FRACTIONAL_SIM_id = find(FRACTIONAL_SIM_arr == 1);
MUR_CONDITION_id = find(MUR_CONDITION_arr == 1);
OPEN_MP_SPACE_id = find(OPEN_MP_SPACE_arr == 1);
TIME_ROW_WISE_id = find(TIME_ROW_WISE_arr == 1);

FRACTIONAL_SIM_id_n = find(FRACTIONAL_SIM_arr == 0);
MUR_CONDITION_id_n = find(MUR_CONDITION_arr == 0);
OPEN_MP_SPACE_id_n = find(OPEN_MP_SPACE_arr == 0);
TIME_ROW_WISE_id_n = find(TIME_ROW_WISE_arr == 0);

%% Ids to filter by - fixed value of alpha, Nz, Nt
NZ = 5000;
NT = 700;

ALPHA = 0.98; 
eps = 1e-8;

NZ_id = find(Nz_arr == NZ);
NT_id = find(Nt_arr == NT);
alpha_id = find((alpha_arr > ALPHA-eps & alpha_arr < ALPHA+eps));

%% Merge desired conditions

% get required flags
FLAGS_id = intersect3(FRACTIONAL_SIM_id, OPEN_MP_SPACE_id, TIME_ROW_WISE_id);
% FLAGS_id = intersect3(FRACTIONAL_SIM_id_n, OPEN_MP_SPACE_id, TIME_ROW_WISE_id);

% get ids for various alpha / Nz / Nt (for constatn Nz and Nt, alpha and Nt, alpha and Nz) 
ALPHA_DEP_id = intersect3(FLAGS_id, NZ_id, NT_id);
NZ_DEP_id = intersect3(FLAGS_id, NT_id, alpha_id);
NT_DEP_id = intersect3(FLAGS_id, NZ_id, alpha_id);



%% Calculate mean and standard deviation, sort as sie-effect

% % Nz dependent
Nz_plt = Nz_arr(NZ_DEP_id);
time_NZ = sim_time_arr(NZ_DEP_id);
[Nz_plt_unique, time_NZ_mean, time_NZ_std] = mean_result(Nz_plt, time_NZ);

% % Nt dependent 
Nt_plt = Nt_arr(NT_DEP_id);
time_NT = sim_time_arr(NT_DEP_id);
[Nt_plt_unique, time_NT_mean, time_NT_std] = mean_result(Nt_plt, time_NT);

% % alpha dependent 
alpha_plt = alpha_arr(ALPHA_DEP_id);
time_ALPHA = sim_time_arr(ALPHA_DEP_id);
[alpha_plt_unique, time_ALPHA_mean, time_ALPHA_std] = mean_result(alpha_plt, time_ALPHA);


plot_Nt_dep = false;
plot_Nz_dep = false;
plot_alp_dep = false;

%% Plot mean value

plot_Nt_dep = true;
% plot_Nz_dep = true;


% Nz dependent   
if plot_Nz_dep
%     scatter(Nz_plt, time_NZ);
%     hold on

    figure();
    scatter(Nz_plt_unique, time_NZ_mean);
    xlabel("Nz")
    ylabel("czas symulacji [s]")
    title("Nt=" + NT + ", \alpha=" + ALPHA);
end

% % Nt dependent 
if plot_Nt_dep
%     scatter(Nt_plt, time_NT);
%     hold on

    figure();
    scatter(Nt_plt_unique, time_NT_mean);
%       semilogy(Nt_plt_unique, time_NT_mean, "o");
%     errorbar(Nt_plt_unique, time_NT_mean, time_NT_std, 'o')

    xlabel("Nt");
    ylabel("czas symulacji [s]");
    title("Nz=" + NZ + ", \alpha=" + ALPHA);
end

% % alpha dependent 
if plot_alp_dep
%     scatter(alpha_plt, time_ALPHA);
%     hold on

    figure()
    scatter(alpha_plt_unique, time_ALPHA_mean);
    xlabel("\alpha")
    ylabel("czas symulacji [s]")
    title("Nz=" + NZ + ", Nt=" + NT);
end

%% intersect3 function
function [Com,ia,ib,ic] = intersect3(A,B,C)
    [C1,ia,ib] = intersect(A,B);
    [Com,ic1,ic] = intersect(C1,C);
    %~ ic is okay
    ia = ia(ic1);
    ib = ib(ic1);
end



%% compute mean (average)
% param var_plt - var dependent values (Nt_plt, Nz_plt, ...)
% param time_var - simulation time for this values
% return:
%   unique_var - array of unique values of var_plt (var_plt without repeating values and sorted)
%   time_var_mean - array of sim times corresponding to unique_var (mean of
%                   different sim times for specific value from unique_var)
%   time_var_std - as above, bu standard deviation
function [unique_var, time_var_mean, time_var_std] = mean_result(var_plt, time_var)
    unique_var = unique(var_plt);  % uniques values of variable - ' to get row vector
    time_var_mean = zeros(size(unique_var));
    time_var_std = zeros(size(unique_var));
    
    for i = 1:length(unique_var) % for every unique value of var_plt ...
        u = unique_var(i);
        u_ids = find(var_plt == u); % ..find indices for this value
        time_var_mean(i) = mean(time_var(u_ids));  % calculate average time, for different sim times for this value of var
        time_var_std(i) = std(time_var(u_ids));
        size(time_var(u_ids), 1) % optionall - print on how many average
    end
end

