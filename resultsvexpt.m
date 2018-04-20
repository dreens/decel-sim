function resultsvexpt(rsf,data)

figure; hold on
scale = 1/60;
colors = [...
    0.2081    0.1663    0.5292;
    0.0582    0.4677    0.8589;
    0.0282    0.6663    0.7574;
    0.4783    0.7489    0.4877;
    0.9264    0.7256    0.2996;
    0 0 0];
speeds = [800 500 200 100 50 37];
for i=1:length(rsf)
    r = rsf(i);
    diamL = 2e-3; 
    distL = 10e-3 + r.pos(1,3) - r.phase*5e-3/180;
    tof = zeros(1,3500);
    times = (1e-7)*(1:length(tof));
    for j=1:length(tof)
        t = times(j);
        zsq = (r.pos(:,3)-distL + t*r.vel(:,3)).^2;
        xsq = (r.pos(:,1)+t*r.vel(:,1)).^2;
        ya = abs(r.pos(:,2)+t*r.vel(:,2));
        tof(j) = sum( zsq + xsq < diamL^2/4 & ya < 1.5e-3);
    end
    c = colors(speeds==r.finalvz,:);
    if length(uniquetol(r.endphases))>1
        l = '-';
    else
        l = '--';
    end
    shift = -(r.finalvz==800)*20;
    plot(times*1e6+shift,tof*scale,'Color',c,'LineStyle',l)
end

for i=1:length(data)
    d = data(i);
    c = colors(d.speed==speeds,:);
    if d.new
        m = '.'; ms = 18;
    else
        m = 'o'; ms = 6;
    end
    shift = -(d.speed==800)*20;
    errorbar(d.XData+shift,d.YData,d.UData,'LineStyle','none','Color',c,'Marker',m,'MarkerSize',ms)
end


end