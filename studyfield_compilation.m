v = VideoWriter('synch_trap_pggg.gif');
v.FrameRate = 1;
open(v);
for i=0:10:90
    h = study2fields('longdecel','singlerod',i);
    drawnow
    %h = studyfield('longdecel',i);
    writeVideo(v,getframe(h));
end
close(v)

