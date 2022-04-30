% Read results

dz = 7.000000e-008;
dt = 2.311599e-017;

filename = "Ex.bin";
fptr=fopen(filename);
spatial_temporal_dimensions =fread(fptr,2,'uint');
Nz=spatial_temporal_dimensions(1);
Nt=spatial_temporal_dimensions(2);
field = fread(fptr, [Nz, Nt],'double');

field_at_t = field(1:end, t);