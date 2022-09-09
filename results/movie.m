function movie(filename)

fptr=fopen(filename);
spatial_temporal_dimensions =fread(fptr,2,'uint');

Nz=spatial_temporal_dimensions(1);
Nt=spatial_temporal_dimensions(2);
dz = fread(fptr,1,'double');
dt = fread(fptr,1,'double');
alpha = fread(fptr,1,'double');
field = fread(fptr, [Nt, Nz],'double');
fclose(fptr);

writerObj = VideoWriter('FieldInTime.avi');
writerObj.FrameRate=10;
open(writerObj);

figure(2);
ylim([min(field(:)) max(field(:))]);

crop = false;
% crop = true;
% Nz_crop = 450;

for t = 1:50:Nt
    figure(2);
    z_arr = (1:Nz).*dz - dz; % -dz to start x-axis from 0
    field_at_t = field(t, 1:end); % field at specific time

    if crop
        z_arr = z_arr(1:Nz_crop);
        field_at_t = field_at_t(1:Nz_crop);
    end
    
    % plot
    plot(z_arr, field_at_t);
    xlabel('z [m]');
    ylabel("Ex [V/m]");
    y_l = get(gca, 'YLim');
    ylim([-max(abs(y_l)), max(abs(y_l))]); % to make zero axis at the center
    xlim([z_arr(1), z_arr(end)])
    title(sprintf('Pole Ex w %e s', (t-1)*dt));  % -1 to start from 0 (el 1 in Matalb - it is initial condition -> 0)
    if t == 1 % scale image if needed
%         sc_h
    end
    % Save frame video
    frame = getframe(gcf);
    writeVideo(writerObj, frame);
    
%     % Save frame to file
%     save_filename = sprintf('t_%d.png', t-1); % -1 to start from 0 (el 1 in Matalb - it is initial condition -> 0)
%     saveas(gcf, save_filename);
    
%     % Save with color
%     save_filename = sprintf('t_%d', t-1); % -1 to start from 0 (el 1 in Matalb - it is initial condition -> 0)
%     saveas(gcf,save_filename, 'epsc');
       
%     % Save in blank an white  
%     save_filename = sprintf('t_%d.eps', t-1); % -1 to start from 0 (el 1 in Matalb - it is initial condition -> 0)
%     saveas(gcf, save_filename);
end
close(writerObj);


end