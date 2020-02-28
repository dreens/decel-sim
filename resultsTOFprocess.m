%% Plots results
function rs = resultsTOFprocess(rs)
    for i=1:length(rs)
        r = rs(i);
        diamL = 1e-3; 
        distL = 10e-3 + 333*5e-3;
        tof = zeros(1,6000);
        times = (1e-7)*(1:length(tof));
        for j=1:length(tof)
            t = times(j);
            zsq = (r.pos(:,3)-distL + t*r.vel(:,3)).^2;
            xsq = (r.pos(:,1)+t*r.vel(:,1)).^2;
            ya = abs(r.pos(:,2)+t*r.vel(:,2));
            tof(j) = sum( zsq + xsq < diamL^2/4 & ya < 1e-3);
        end
        rs(i).tofpeak = max(tof);
        rs(i).tofarea = trapz(times,tof);
    end
end