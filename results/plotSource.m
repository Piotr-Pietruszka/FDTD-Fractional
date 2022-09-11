function [t_sec, source_in_time]=plotSource(filename)


% Load source from file
fptr=fopen(filename);
spatial_temporal_dimensions =fread(fptr,2,'uint');

Nz = spatial_temporal_dimensions(1);
Nt = spatial_temporal_dimensions(2);
dz = fread(fptr,1,'double');
dt = fread(fptr,1,'double');
alpha = fread(fptr,1,'double');
k_bound = fread(fptr,1,'int');
source_in_time = fread(fptr, [Nz, Nt],'double');
fclose(fptr);


% Alternatively - chosen point at domain
% ------------
% fptr=fopen("Ex.bin");
% spatial_temporal_dimensions =fread(fptr,2,'uint');
% Nz=spatial_temporal_dimensions(1);
% Nt=spatial_temporal_dimensions(2);
% dz = fread(fptr,1,'double');
% dt = fread(fptr,1,'double');
% alpha = fread(fptr,1,'double');
% k_bound = fread(fptr,1,'int');
% field = fread(fptr, [Nt, Nz],'double');
% source_in_time = field(1:end, int32(0.1e-6/dz)+1);
% % source_in_time = field(1:end, int32(0.1e-6/dz)+int32(0e-6/dz)++1);
% fclose(fptr);
% ------------

% crop source in time
% Nt=1500;
% source_in_time = source_in_time(1:Nt);
% max(source_in_time)


t_sec=(0:1:Nt-1)*dt; % t in secs
% t_sec=(1:Nt).*dt - dt;; % t in secs
t=(0:1:Nt-1); % cells
% size(source_in_time)
% size(t)

figure(1)
plot(t_sec, source_in_time); 
title('Pole �r�d�owe');
xlabel('t (s)');
ylabel("Ex [V/m]");
% xlim([0, 2.5e-14]); % "crop" - zxlim

% Spectrum of the source
y = fft(source_in_time, Nt)*2/Nt;
y = fftshift(y); % shift f transform
fs = 1/dt;
fshift = (-Nt/2:Nt/2-1)*(fs/Nt); % domain in Hz
size(fshift)

% crop frequency - get only right half
% ix1 = find(fshift>=0, 1)
% fshift = fshift(ix1:end);
% y = y(ix1:end);

% % specific range
% ix1 = find(fshift>=3e14, 1)
% ix2 = find(fshift>=9e14, 1)
% fshift = fshift(ix1:ix2);
% y = y(ix1:ix2);

figure(2)
plot(fshift, abs(y))
title('Spektrum �r�d�a');
xlabel('f (Hz)');


assignin('base','t_sec',t_sec)
assignin('base','source_in_time',source_in_time)
end