v = VideoWriter('synch_trap2');
v.FrameRate = 1;
open(v);
for p=0:5:90
    h = study2fields('longdecel','longdecel',p-180,-p,p);
    drawnow
    %h = studyfield('longdecel',i);
    writeVideo(v,getframe(h));
end
close(v)
close all
