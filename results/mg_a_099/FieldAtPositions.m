function [f2, t_arr]=FieldAtPositions(filename, z_points)

fptr=fopen(filename);

   
spatial_temporal_dimensions = fread(fptr,2,'uint');
Nz = spatial_temporal_dimensions(1);
Nt = spatial_temporal_dimensions(2);
dz = fread(fptr,1,'double');
dt = fread(fptr,1,'double');
alpha = fread(fptr,1,'double');
field = fread(fptr, [Nt, Nz],'double');
fclose(fptr);

% return all fields for given positions:
% z1   z2   z3
% _    _    _  t=0
% _    _    _  t=1
% _    _    _  t=2
% .    .    .   .
f2 = zeros(Nt, length(z_points)); % Nt - no of rows, no fields = no. columns 

t_arr = (1:Nt).*dt - dt; % -dt to start x-axis from 0

for i = 1:length(z_points)
    z = z_points(i) + 1; % +1: give position as in C (from 0) -> position in Matlab (from 1)
    field_at_z = field(1:end, z);  % field at specific time
        
    % plot
%     plot(t_arr, field_at_z);
%     xlabel('t [s]');
%     ylabel("Ex [V/m]");
%     xlim([t_arr(1), t_arr(end)]);
%     title(sprintf('Pole Ex w punkcie %e m', double(z-1)*dz)); % -1 to start from 0 (el 1 in Matalb - it is boundary -> 0)
    
%     save_filename = sprintf('z_%d.png', z-1); % -1 to start from 0 (el 1 in Matalb - it is boundary -> 0)
%     saveas(gcf, save_filename);
    
%     assignin('base','t_arr',t_arr);
%     assignin('base','field_at_z',field_at_z);
    f2(:, i) = field_at_z; % column - get field for specific position
end


end
