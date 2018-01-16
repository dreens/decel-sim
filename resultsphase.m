%% Makes some phase space plots
function resultsphase(rs)
    figure;
    run = rs(1);
    lLZ = setdiff(run.whichmolload,run.whichmol);
    lLoad = setdiff(1:10000,run.whichmolload);
    nl = run.whichmol;
    x =  run.posi(nl,1)*1e3;
    y =  run.posi(nl,2)*1e3;
    z =  run.posi(nl,3)*1e3;
    vx = run.veli(nl,1);
    vy = run.veli(nl,2);
    vz = run.veli(nl,3);
    r = sqrt(x.^2+y.^2);
    t = atan(y./x);
    vr = (x.*vx+y.*vy)./r;
    vt = (x.*vy-y.*vx)./r;
    xd =  run.posi(lLoad,1)*1000;
    yd =  run.posi(lLoad,2)*1000;
    zd =  run.posi(lLoad,3)*1000;
    vxd = run.veli(lLoad,1);
    vyd = run.veli(lLoad,2);
    vzd = run.veli(lLoad,3);
    rd = sqrt(xd.^2+yd.^2);
    td = atan(yd./xd);
    vrd = (xd.*vxd+yd.*vyd)./rd;
    vtd = (xd.*vyd-yd.*vxd)./rd;
    xl =  run.posi(lLZ,1)*1e3;
    yl =  run.posi(lLZ,2)*1e3;
    zl =  run.posi(lLZ,3)*1e3;
    vxl = run.veli(lLZ,1);
    vyl = run.veli(lLZ,2);
    vzl = run.veli(lLZ,3);
    rl = sqrt(xl.^2+yl.^2);
    tl = atan(yl./xl);
    vrl = (xl.*vxl+yl.*vyl)./rl;
    vtl = (xl.*vyl-yl.*vxl)./rl;
    m = 3;
    subplot(2,3,1)
    plot(xd,yd,'b.','MarkerSize',m); hold on;
    plot(x,y,'g.','MarkerSize',m);
    plot(xl,yl,'r.','MarkerSize',m);
    xlabel('x (mm)')
    ylabel('y (mm)')
    axis equal
    title('x-y')
    subplot(2,3,2)
    plot(rd,zd,'b.','MarkerSize',m); hold on;
    plot(r,z,'g.','MarkerSize',m);
    plot(rl,zl,'r.','MarkerSize',m);
    xlabel('r (mm)')
    ylabel('z (mm)')
    axis equal
    title('r-z')
    title({'Initial Cylindrical Phase Space, Labeled by Outcome' 'r-z'})
    subplot(2,3,3)
    plot(td,zd,'b.','MarkerSize',m); hold on;
    plot(t,z,'g.','MarkerSize',m);
    plot(tl,zl,'r.','MarkerSize',m);
    xlabel('\theta (rad)')
    ylabel('z (mm)')
    axis equal
    title('\theta-z')
    subplot(2,3,4)
    plot(vrd,vtd,'b.','MarkerSize',m); hold on;
    plot(vr,vt,'g.','MarkerSize',m);
    plot(vrl,vtl,'r.','MarkerSize',m);
    xlabel('v_r (m/s)')
    ylabel('v_\theta (m/s)')
    axis equal
    title('v_r-v_\theta')
    subplot(2,3,5)
    plot(vrd,vzd,'b.','MarkerSize',m);hold on;
    plot(vr,vz,'g.','MarkerSize',m); 
    plot(vrl,vzl,'r.','MarkerSize',m);
    xlabel('v_r (m/s)')
    ylabel('v_z (m/s)')
    axis equal
    title({'Blue = Loading Loss    ' 'Green = Trapped          ' 'Red = Majoranna Loss' '' 'v_r-v_z'})
    subplot(2,3,6)
    plot(vtd,vzd,'b.','MarkerSize',m); hold on;
    plot(vt,vz,'g.','MarkerSize',m);
    plot(vtl,vzl,'r.','MarkerSize',m);
    axis equal
    xlabel('v_\theta (m/s)')
    ylabel('v_z (m/s)')
    title('v_\theta-v_z')
    
    figure;
    plot(rd,vtd,'b.','MarkerSize',m); hold on;
    plot(r,vt,'g.','MarkerSize',m);
    plot(rl,vtl,'r.','MarkerSize',m);
    xlabel('r (mm)')
    ylabel('v_\theta (m/s)')
end
