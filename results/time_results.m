
%% Read data from file
filename = "results.txt";
% filename = "backup\\results_Nt_all.txt";
% filename = "backup\\results_Nz_all.txt";
% filename = "backup\\results_alpha_all.txt";
% filename = "backup\\results_column_Nt_all.txt";


d=importdata(filename);




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

%% Ids to filter by - flags and fixed values of alpha, Nz, Nt

% get required flags
FLAGS_id = intersect3(FRACTIONAL_SIM_id, OPEN_MP_SPACE_id, TIME_ROW_WISE_id); TIME_ROW_WISE = 1; FRACTIONAL_SIM = 1;
% FLAGS_id = intersect3(FRACTIONAL_SIM_id, OPEN_MP_SPACE_id, TIME_ROW_WISE_id_n); TIME_ROW_WISE = 0; FRACTIONAL_SIM = 1;

% FLAGS_id = intersect3(FRACTIONAL_SIM_id_n, OPEN_MP_SPACE_id, TIME_ROW_WISE_id); FRACTIONAL_SIM = 0;

NZ = 7000;
NT = 3000;
T = 1e-13;
ALPHA = 0.98; 
if ~FRACTIONAL_SIM
    ALPHA = 1; 
end

eps = 1e-8;

NZ_id = find(Nz_arr == NZ); % Nz ids for specific value
NT_id = find(Nt_arr == NT); % Nt ids for specific value
alpha_id = find((alpha_arr > ALPHA-eps & alpha_arr < ALPHA+eps)); % alpha ids for alpha from interval (alpha - float)
T_id = find(abs(T_arr - T) < eps); % T ids fo T from interval (T - float)


%% Merge desired conditions

% get ids for various alpha / Nz / Nt (for constatn Nz and Nt, alpha and Nt, alpha and Nz) 
ALPHA_DEP_id = intersect3(FLAGS_id, NZ_id, T_id);
NZ_DEP_id = intersect3(FLAGS_id, NT_id, alpha_id);
NT_DEP_id = intersect3(FLAGS_id, NZ_id, alpha_id);


plot_Nt_dep = false;
plot_Nz_dep = false;
plot_alp_dep = false;


%% Calculate mean and standard deviation, sort as sie-effect

% plot_Nz_dep = true;
plot_Nt_dep = true;
% plot_alp_dep = true;

% % Nz dependent
if plot_Nz_dep
    Nz_plt = Nz_arr(NZ_DEP_id);
    time_NZ = sim_time_arr(NZ_DEP_id);
    [Nz_plt_unique, time_NZ_mean, time_NZ_std] = mean_result(Nz_plt, time_NZ);
end

% % Nt dependent 
if plot_Nt_dep
    Nt_plt = Nt_arr(NT_DEP_id);
    time_NT = sim_time_arr(NT_DEP_id);
    [Nt_plt_unique, time_NT_mean, time_NT_std] = mean_result(Nt_plt, time_NT);
end

% % alpha dependent 
if plot_alp_dep
    alpha_plt = alpha_arr(ALPHA_DEP_id);
    time_ALPHA = sim_time_arr(ALPHA_DEP_id);
    [alpha_plt_unique, time_ALPHA_mean, time_ALPHA_std] = mean_result(alpha_plt, time_ALPHA);
end



%% Plot mean value




% Nz dependent   
if plot_Nz_dep
%     scatter(Nz_plt, time_NZ);
%     hold on

%     figure();
    scatter(Nz_plt_unique, time_NZ_mean);
    xlabel("Nz")
    ylabel("czas wykonywania symulacji [s]");
    title("Czas wykonywania, \alpha=" + ALPHA);
    
    for i = 1:size(Nz_plt_unique, 1)
        print_latex(time_NZ_mean(i), time_NZ_std(i), NT, Nz_plt_unique(i), ALPHA, TIME_ROW_WISE);
    end
end

% % Nt dependent 
if plot_Nt_dep
%     scatter(Nt_plt, time_NT);
%     hold on

%     figure();
    scatter(Nt_plt_unique, time_NT_mean);
%       semilogy(Nt_plt_unique, time_NT_mean, "o");
%     errorbar(Nt_plt_unique, time_NT_mean, time_NT_std, 'o')

    xlabel("Nt");
    ylabel("czas wykonywania symulacji [s]");
%     title("Nz=" + NZ + ", \alpha=" + ALPHA);
    title("Czas wykonywania, \alpha=" + ALPHA);
    
    for i = 1:size(Nt_plt_unique, 1)
        print_latex(time_NT_mean(i), time_NT_std(i), Nt_plt_unique(i), NZ, ALPHA, TIME_ROW_WISE);
    end
end

% % alpha dependent 
if plot_alp_dep
%     scatter(alpha_plt, time_ALPHA);
%     hold on

%     figure();
    scatter(alpha_plt_unique, time_ALPHA_mean);
    xlabel("\alpha")
    ylabel("czas wykonywania symulacji [s]");
    title("Nz=" + NZ + ", T=" + T + "[s]");

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

%% print in latex array foramt
% all param - scalars
% param time_mean - mean time (in seconds)
% param time_std - standard deviation (in seconds)
% param Nt - no time steps
% param Nz - no sptaial cells
% param alpha - value of alpha 
% param TIME_ROW_WISE - 1 - row wize, 0 - columnwise
function v = print_latex(time_mean, time_std, Nt, Nz, alpha, TIME_ROW_WISE)
    time_row_column_array = ["kolumnowe", "wierszowe"];
%      4000  &  1000  &   0.98   &  wierszowe  &  112.33333   &   1  \% \\ \hline

    l = Nz + "  &  " + Nt + "  &   " + alpha + "   &  " + time_row_column_array(TIME_ROW_WISE+1) + ...
        "  &  " + round(time_mean, 2) + "   &   " + round( time_std / time_mean*100, 2) + "\% \\ \hline";
    disp(l);
    
end

