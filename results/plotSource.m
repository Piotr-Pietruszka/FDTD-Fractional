function plotSource(filename)


% dt = 3.302285e-017;
dt = 2.311599e-017;
% Load source from file
fptr=fopen(filename);
spatial_temporal_dimensions =fread(fptr,2,'uint');

Nz=spatial_temporal_dimensions(1);
Nt=spatial_temporal_dimensions(2);
source_in_time = fread(fptr, [Nz, Nt],'double');


% crop source in time
Nt=700;
source_in_time = source_in_time(1:Nt);


t=(0:1:Nt-1)*dt; % t in secs
% size(source_in_time)
% size(t)

plot(t, source_in_time); 
title('E field amplitude');
xlabel('t (s)');

% Spectrum of the source
figure(2)
y = fft(source_in_time, Nt);
y = fftshift(y); % shift f transform
fs = 1/dt;
fshift = (-Nt/2:Nt/2-1)*(fs/Nt); % domain in Hz
% f = (0:length(y)-1)*fs/length(y);
% plot(f, abs(y));
size(fshift)
% crop frequency - get second only right half
ix1 = find(fshift>=3e14, 1)
ix2 = find(fshift>=9e14, 1)
fshift = fshift(ix1:ix2);
y = y(ix1:ix2);

plot(fshift, abs(y)/max(abs(y))*0.13)

% plot(source_in_time);
% title('E field amplitude');
% xlabel('t (steps)');


end