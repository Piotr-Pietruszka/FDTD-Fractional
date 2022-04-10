function movie(filename)

fptr=fopen(filename);
spatial_temporal_dimensions =fread(fptr,2,'uint');

Nz=spatial_temporal_dimensions(1);
Nt=spatial_temporal_dimensions(2);
field = fread(fptr, [Nz, Nt],'double');


writerObj = VideoWriter('YourAVI.avi');
writerObj.FrameRate=1;
open(writerObj);
for t = 1:400:Nt
    field_at_t = field(1:end, t);
    plot(field_at_t);
    title('E field amplitude');
    xlabel('z (cells)');
    frame = getframe(gcf);
    writeVideo(writerObj, frame);
%     plot(field_at_t)
end
close(writerObj);


end