v = VideoWriter('synch_trap_ppmmonly');
v.FrameRate = 1;
open(v);
for p=0:5:90
    h = study2fields('longdecel','ppmm_2mm',-135+p/2,-45+p/2,p);
    drawnow
    %h = studyfield('longdecel',i);
    writeVideo(v,getframe(h));
end
close(v)
close all
