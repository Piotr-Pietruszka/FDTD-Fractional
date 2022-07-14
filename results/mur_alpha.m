
MU_0 = 1.2566371e-6;
EPS_0 = 8.85418781762e-12;
C_CONST = 299792458.0;

% z_m = 6e-6;  % z from left in meters
z_m = (20e-6) / 3

%% PEC
alpha_arr = [];
max_E_arr = [];
folder = "";
folder = "C2/"

for i = 0:20
    %% Read file
    filename = folder + "Ex_PEC_" + i + ".bin";
    fptr=fopen(filename);

    spatial_temporal_dimensions = fread(fptr,2,'uint');
    Nz = spatial_temporal_dimensions(1);
    Nt = spatial_temporal_dimensions(2);
    dz = fread(fptr,1,'double');
    dt = fread(fptr,1,'double');
    alpha = fread(fptr,1,'double');
    % Get Courant coefficient
    fractional_dt_analytical = 2^(1-1/alpha) * (sqrt(EPS_0*MU_0) * dz)^(1/alpha);
    S = dt / fractional_dt_analytical;
    
    field = fread(fptr, [Nt, Nz],'double');
    fclose(fptr);
    z = int32(z_m / dz);
    
    field_at_z = field(1:end, z);
    
    %% Max reflected value for alpha
    alpha_arr(i+1) = alpha;
    
    % czas uzyskac poprze wziecie takiego punktu, ktory jest w odleglosci
    % 1/3 od granicy i 2/3 od zrodla (w z) jako z. Tam wziac pole f
    % W polu maksimum bedzie pik fali padajacej. To moment przyjscia fali
    % padajacej.
    % Po tym samym czasie powinna przyjsc fala odbita (2*1/3 = 2/3)
    % Wiec jako punkt graniczny wybrac ten pomiedzy - wiêc graniczny punkt w czasie =
    % t_pik_1 + t_pik_1/2   
    
    [max_val, t_max] = max( abs(field_at_z) );
%     max_E_arr(i+1) = max( abs(field_at_z(int32(5e-14/dt):end)) ) /  max( abs(field_at_z) );  
    t_slice = int32(1.5*t_max);
    max_E_arr(i+1) = max( abs(field_at_z(t_slice:end)) );  

end
S

%% Plot results
% plot(alpha_arr, max_E_arr, "o")
% grid on;
% xlabel("\alpha")
% ylabel("Odbicie")


%% MUR

alpha_arr_mur = [];
max_E_arr_mur = [];
for i = 0:20
    %% Read file
    filename = folder + "Ex_MUR_" + i + ".bin";
    fptr=fopen(filename);

    spatial_temporal_dimensions = fread(fptr,2,'uint');
    Nz = spatial_temporal_dimensions(1);
    Nt = spatial_temporal_dimensions(2);
    dz = fread(fptr,1,'double');
    dt = fread(fptr,1,'double');
    alpha = fread(fptr,1,'double');
    
    field = fread(fptr, [Nt, Nz],'double');
    fclose(fptr);
    z = int32(z_m / dz);
    
    field_at_z = field(1:end, z);
    
    %% Max reflected value for alpha
    alpha_arr_mur(i+1) = alpha;
    
    [max_val, t_max] = max( abs(field_at_z) );
    t_slice = int32(1.5*t_max);
    max_E_arr_mur(i+1) = max( abs(field_at_z(t_slice:end)) );  

end


%% Plot results - reflection
plot(alpha_arr, max_E_arr_mur ./ max_E_arr , "o")
% plot(alpha_arr, max_E_arr_mur ./ max_E_arr )

grid on;
xlabel("\alpha")
ylabel("Odbicie")
title("Odbicie dla warunku Mura nieca³kowitego rzêdu")
