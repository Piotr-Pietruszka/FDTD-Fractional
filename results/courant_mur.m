

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