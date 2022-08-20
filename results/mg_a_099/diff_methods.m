a=0.99;
b=0.99;

filename = 'Ex_dz_0005.bin';
dz = 0.005e-6;

L_arr = (1:1:130) * 1e-6;

max_fdtd_arr = zeros(size(L_arr));
max_f_arr = zeros(size(L_arr));
max_diff_arr = zeros(size(L_arr));
max_diff_id_arr = zeros(size(L_arr));

max_fdtd_arr_for_diff = zeros(size(L_arr));
max_f_arr_for_diff = zeros(size(L_arr));

i = 1;

all_at_once = true; % whether to read field at all positions at once
% all_at_once = false;
if all_at_once
    [f2, t_arr] = FieldAtPositions(filename, int32(0.1e-6/dz)+int32(L_arr/dz));
end 

for L=L_arr
%    max_f_arr(i) = max(abs(output_signal4));
%    max_fdtd_arr(i) = max(abs(plotFieldAtPositions('Ex.bin', [int32(0.1e-6/dz)+int32(L/dz)])));
   
   % Get arrays in time, for position for 2 methods
   if ~all_at_once
       fdtd_arr = FieldAtPositions(filename, [int32(0.1e-6/dz)+int32(L/dz)]);
       
   else
       fdtd_arr = f2(:, i);
   end
   [ t4, f_arr ] = superluminal_1( a, b, L );

   % get max difference
%    [max_diff_arr(i), max_diff_id_arr(i)] = max(abs(fdtd_arr' - f_arr));
   % get values from arrays for same id as max difference 
%    max_fdtd_arr_for_diff(i) = fdtd_arr(max_diff_id_arr(i));
%    max_fdtd_arr_for_diff(i) = f_arr(max_diff_id_arr(i));

   % Get max value for whole arrays
   max_fdtd_arr(i) = max(abs(fdtd_arr));
   max_f_arr(i) = max(abs(f_arr));
   
   % Dimensions are much different (different T, dt), so this is taken as
   % max diff
   max_diff_arr(i) = abs(max_fdtd_arr(i) - max_f_arr(i));
   i = i+1;
   if mod(i, 20) == 0
       i
   end
end

plot(L_arr+0.1e-6, max_diff_arr./max_f_arr); % +0.1e-6 - source position
xlabel('z [m]');
ylabel("Maksymalny b³¹d");
title("Ró¿nica wyników miêdzy metodami");
% title(sprintf('Pole Ex w punkcie %e m', double(z-1)*dz)); 
    
