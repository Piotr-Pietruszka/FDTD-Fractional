function movie(filename)

fptr=fopen(filename);
spatial_temporal_dimensions =fread(fptr,2,'uint');

Nz=spatial_temporal_dimensions(1);
Nt=spatial_temporal_dimensions(2);
dz = fread(fptr,1,'double');
dt = fread(fptr,1,'double');
field = fread(fptr, [Nt, Nz],'double');
fclose(fptr);

writerObj = VideoWriter('FieldInTime.avi');
writerObj.FrameRate=10;
open(writerObj);

figure(2);
ylim([min(field(:)) max(field(:))]);

for t = 1:100:Nt
% for t = 1:10:600
    figure(2);
    if exist('dz') > 0
        z_arr = (1:Nz).*dz;
        xl = 'z [m]';
    else
        z_arr = (1:Nz);
        xl = 'z [cells]';
    end
    field_at_t = field(t, 1:end);
    plot(z_arr, field_at_t);
    xlabel(xl);
    y_l = get(gca, 'YLim');
    ylim([-max(abs(y_l)), max(abs(y_l))]);
    
    % Title with curretn time (sec / steps) depending on existence of dt
    if exist('dt') > 0
        title(sprintf('E field amplitude at %e s', t*dt));
    else
        title(sprintf('E field amplitude at step %d', t));
    end
    
    % Save frame video
    frame = getframe(gcf);
    writeVideo(writerObj, frame);
    
    % save frame to file
    save_filename = sprintf('t_%d.png', t);
    saveas(gcf, save_filename);
end
close(writerObj);


end