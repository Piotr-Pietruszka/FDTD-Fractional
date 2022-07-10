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

for i = 1:length(z_points)
    z = z_points(i);
    field_at_z = field(1:end, z);
    
    if filename == "Ex.bin" % TEMP
%         field_at_z = field_at_z*0.6599/0.305;
%         a = 1
    end
    % X-axis units - depnd on existence of dt
    if exist('dt') > 0
        t_arr = (1:Nt).*dt;
        xl = 't [s]';
    else
        t_arr = (1:Nt);
        xl = 't [steps]';
    end
    plot(t_arr, field_at_z);
    xlabel(xl);
    
    % Title depends on existence of dz
    if exist('dz') > 0
        title(sprintf('E field amplitude at %e m', double(z)*dz));
    else
        title(sprintf('E field amplitude at cell %d', z));
    end
    
    
    save_filename = sprintf('z_%d.png', z);
    saveas(gcf, save_filename);
end


end
