function field=readField(filename)


fptr=fopen(filename);
spatial_temporal_dimensions =fread(fptr,2,'uint');

Nz=spatial_temporal_dimensions(1);
Nt=spatial_temporal_dimensions(2);
    
field = fread(fptr, [Nz, Nt],'double');


mesh(field);
view(0,90);
caxis([-1 1]);
size(field)


title('E field amplitude');
xlabel('t (steps)');
ylabel('z (cells)');


% % % Linear
colorbar;
colormap(hsv);


fclose(fptr);
