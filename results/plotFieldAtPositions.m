function plotFieldAtPositions(filename, z_points)

fptr=fopen(filename);
spatial_temporal_dimensions =fread(fptr,2,'uint');

Nz=spatial_temporal_dimensions(1);
Nt=spatial_temporal_dimensions(2);
field = fread(fptr, [Nt, Nz],'double');

dt=3.302285e-017;
for i = 1:length(z_points)
    z = z_points(i);
    field_at_z = field(1:end, z);
    plot((1:Nt).*dt, field_at_z);
    title(strcat('E field amplitude at cell: ', num2str(z)));
    xlabel('t (steps)');
    
    save_filename = strcat(num2str(z), "_z.png");
    saveas(gcf, save_filename);
end


end
