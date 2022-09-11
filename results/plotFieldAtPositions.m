function field_at_z=plotFieldAtPositions(filename, z_points)

fptr=fopen(filename);

   
spatial_temporal_dimensions = fread(fptr,2,'uint');
Nz = spatial_temporal_dimensions(1);
Nt = spatial_temporal_dimensions(2);
dz = fread(fptr,1,'double');
dt = fread(fptr,1,'double');
alpha = fread(fptr,1,'double');
k_bound = fread(fptr,1,'int');
field = fread(fptr, [Nt, Nz],'double');
fclose(fptr);

% crop = true;
crop = false;

for i = 1:length(z_points)
    z = z_points(i) + 1; % +1: give position as in C (from 0) -> position in Matlab (from 1)
    field_at_z = field(1:end, z);  % field at specific time
    
    t_arr = (1:Nt).*dt - dt; % -dt to start x-axis from 0

    if crop
        field_at_z = field_at_z(150:end);
        t_arr = t_arr(150:end);
    end
    
    % plot
    figure()
    plot(t_arr, field_at_z);
    xlabel('t [s]');
    ylabel("Ex [V/m]");
    xlim([t_arr(1), t_arr(end)]);
    title(sprintf('Pole Ex w punkcie %e m', double(z-1)*dz)); % -1 to start from 0 (el 1 in Matalb - it is boundary -> 0)
    
%     save_filename = sprintf('z_%d.png', z-1); % -1 to start from 0 (el 1 in Matalb - it is boundary -> 0)
%     saveas(gcf, save_filename);
    
    size(field_at_z);
    min(field_at_z);
    
%     assignin('base','t_arr',t_arr);
%     assignin('base','field_at_z',field_at_z);
end


end
