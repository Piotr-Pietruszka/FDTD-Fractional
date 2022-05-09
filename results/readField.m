function field=readField(filename)


fptr=fopen(filename);
spatial_temporal_dimensions =fread(fptr,2,'uint');

Nz=spatial_temporal_dimensions(1);
Nt=spatial_temporal_dimensions(2);
dz = fread(fptr,1,'double');
dt = fread(fptr,1,'double');    
field = fread(fptr, [Nt, Nz],'double');
% field = field';

mesh(field);
view(0,90);
caxis([-1 1]);
size(field)


title('E field amplitude');
xlabel('z (cells)');
ylabel('t (steps)');

% % % Linear
colorbar;
colormap(hsv);


fclose(fptr);
