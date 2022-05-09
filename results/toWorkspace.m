% Read results


filename = "Ex.bin";
fptr=fopen(filename);
spatial_temporal_dimensions =fread(fptr,2,'uint');
Nz=spatial_temporal_dimensions(1);
Nt=spatial_temporal_dimensions(2);
dz = fread(fptr,1,'double');
dt = fread(fptr,1,'double');
field = fread(fptr, [Nt, Nz],'double');
fclose(fptr);

field_at_t = field(t, 1:end);

