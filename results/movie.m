function movie(filename)

fptr=fopen(filename);
spatial_temporal_dimensions =fread(fptr,2,'uint');

Nz=spatial_temporal_dimensions(1);
Nt=spatial_temporal_dimensions(2);
field = fread(fptr, [Nz, Nt],'double');


writerObj = VideoWriter('FieldInTime.avi');
writerObj.FrameRate=10;
open(writerObj);

figure(2);
ylim([min(field(:)) max(field(:))]);

for t = 1:100:Nt
    field_at_t = field(1:end, t);
    plot(field_at_t);
%     ylim([min(field(:)) max(field(:))]);
    y_l = get(gca, 'YLim');
    ylim([-max(abs(y_l)), max(abs(y_l))]);
    title('E field amplitude');
    xlabel('z (cells)');
    frame = getframe(gcf);
    writeVideo(writerObj, frame);
%     plot(field_at_t)
end
close(writerObj);


end