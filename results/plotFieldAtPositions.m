function field_at_z=plotFieldAtPositions(filename, z_points)

fptr=fopen(filename);


    
spatial_temporal_dimensions = fread(fptr,2,'uint');
Nz = spatial_temporal_dimensions(1);
Nt = spatial_temporal_dimensions(2);
dz = fread(fptr,1,'double');
dt = fread(fptr,1,'double');
alpha = fread(fptr,1,'double');
field = fread(fptr, [Nt, Nz],'double');
fclose(fptr);

% crop = true;
crop = false;

for i = 1:length(z_points)
    z = z_points(i);
    field_at_z = field(1:end, z);
    

    % X-axis units - depnd on existence of dt
    if exist('dt') > 0
        t_arr = (1:Nt).*dt - dt;
        xl = 't [s]';
    else
        t_arr = (1:Nt) - 1;
        xl = 't [steps]';
    end
    
    if crop
        field_at_z = field_at_z(150:end);
        t_arr = t_arr(150:end);
    end
    plot(t_arr, field_at_z);
    xlabel(xl);
    ylabel("Ex [V/m]");
    
    % Title depends on existence of dz
    if exist('dz') > 0
        title(sprintf('Pole Ex w punkcie %e m', double(z-1)*dz));
    else
        title(sprintf('Pole Ex w komórce %d', z));
    end
    
    
    save_filename = sprintf('z_%d.png', z);
    saveas(gcf, save_filename);
    
    size(field_at_z)
    min(field_at_z)
end


end
